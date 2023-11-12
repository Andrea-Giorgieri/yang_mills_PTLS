#ifndef PARALLEL_TEMPERING_C
#define PARALLEL_TEMPERING_C

#include"../include/macro.h"

#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/random.h"
#include"../include/function_pointers.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"

// swaps are parallelized, evaluation of swap probabilities is parallelized
void swap(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
				 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)
	{
	int aux_i, i, j, err, num_swaps, is_even, is_even_first;
	long k, s, num_even, num_odd, num_swaps_1, num_swaps_2;
	double *metro_swap_prob;
	
	// N_replica_pt - 1 is the total number of swaps
	num_swaps = ((param->d_N_replica_pt)-1);
	
	// define auxiliary array to store metropolis swap probabilities
	err=posix_memalign( (void **) &(metro_swap_prob), (size_t) DOUBLE_ALIGN, (size_t) num_swaps * sizeof(double));
	if(err!=0)
		{
		fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	
	// set all probabilities to 0
	for(k=0; k<num_swaps; k++)
		metro_swap_prob[k]=0.0;

	is_even = num_swaps % 2;                   // check if num_swaps is even or not
	num_even = (long) ((num_swaps+is_even)/2); // number or swaps for even replica
	num_odd  = (long) ((num_swaps-is_even)/2); // number of swaps for odd replica
	
	// to be sure detailed balance is satisfied, choose randomly whether to swap first odd or even copies
	
	if( casuale() < 0.5 ) // first swap all even copies, then all odd copies 
		{
		is_even_first=0;
		num_swaps_1=num_even;
		num_swaps_2=num_odd;
		}
	else // first swap all odd copies, then all even copies 
		{
		is_even_first=1;
		num_swaps_1=num_odd;
		num_swaps_2=num_even;
		}

	// first group of swaps	
	
	// compute action differences
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:metro_swap_prob[:num_swaps]) private(s,aux_i,i,j)
	#endif
	for(s=0;s<((num_swaps_1)*(swap_rectangle->d_vol_rect));s++)
		{
		long n = s%(swap_rectangle->d_vol_rect);
		long r = swap_rectangle->rect_sites[n]; // action changes only on the highest level of the defect
		int l = (int) ( (s-n) / (swap_rectangle->d_vol_rect) );
		int a = 2*l+is_even_first; // labels of replica
		int b = a+1;
		
		// swaps are done for all couples (i,j) where j !=i => six couples
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
				{
				// contribution to action difference between replicas a and b of site r on plane (i,j)
				metro_swap_prob[a] += delta_action_swap(GC, geo, param, r, i, j, a, b);
				}		
		}
	
	// do the swaps
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(k)
	#endif
	for(k=0;k<num_swaps_1; k++)
		{
		int a=2*((int)k)+is_even_first;
		int b=a+1;
		metro_swap_prob[a] = exp(-metro_swap_prob[a]); // metropolis swap probability = exp{ - (swapped action - unswapped action) }
		metropolis_single_swap(GC, a, b, metro_swap_prob[a], acc_counters);
		}

	// second group of swaps
	
	is_even_first=1-is_even_first; // used to pass from swapping even copies to odd copies and viceversa
	
	// compute action differences
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) reduction(+:metro_swap_prob[:num_swaps]) private(s,aux_i,i,j)
	#endif
	for(s=0;s<((num_swaps_2)*(swap_rectangle->d_vol_rect));s++)
		{
		long n = s%(swap_rectangle->d_vol_rect);
		long r = swap_rectangle->rect_sites[n]; // action changes only in the first neighborhood of the defect (having swapped also twist factors) 
		int l = (int) ( (s-n) / (swap_rectangle->d_vol_rect) );
		int a = 2*l+is_even_first; // labels of replica
		int b = a+1;
		
		// swaps are done for all couples (i,j) where j !=i => six couples
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
				{
				// contribution to action difference between replicas a and b of site r on plane (i,j)
				metro_swap_prob[a] += delta_action_swap(GC, geo, param, r, i, j, a, b);
				}
		}
	
	// do the swaps
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(k)
	#endif
	for(k=0;k<num_swaps_2; k++)
		{
		int a=2*((int)k)+is_even_first;
		int b=a+1;
		metro_swap_prob[a] = exp(-metro_swap_prob[a]); // metropolis swap probability = exp{ - (swapped action - unswapped action) }
		metropolis_single_swap(GC,a,b,metro_swap_prob[a],acc_counters); // metropolis step
		}
	
	// free aux array
	free(metro_swap_prob);
	}

double delta_action_swap(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param,
                         long const r, int const i, int const j, int const a, int const b)
	{
	double re_tr_plaq_a, re_tr_plaq_b, delta;

	// plaquettes, including twist factors in function plaquettep
	re_tr_plaq_a = plaquettep(&(GC[a]),geo,param,r,i,j); // (Re Tr plaq_a(r,i,j) )/N_c , replica label = a
	re_tr_plaq_b = plaquettep(&(GC[b]),geo,param,r,i,j); // (Re Tr plaq_b(r,i,j) )/N_c , replica label = b

	// (swapped action - unswapped action) = delta_beta * delta_plaq (twist factors swapped, otherwise delta_beta -> beta_a*Z_a-beta_b*Z_b )
	delta = (GC[a].B[r] - GC[b].B[r]) * (re_tr_plaq_a - re_tr_plaq_b);
	
	return delta;
	}

// swaps are serial, evaluation of swap probability is parallelized (use this version of 'swap' if gcc_version < 6.0 or icc_version < 14.0)
/*
void swap(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
				 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters)
  {
	int aux_i, i, j, num_swaps, a, b;
	long n;
	double metro_swap_prob = 0.0;
	
	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };
	
	// N_replica_pt - 1 is the total number of swaps
	num_swaps = ((param->d_N_replica_pt)-1);

	// swaps are done for all couples (i,j) where i=defect_dir and j !=i => three couples
	i=param->d_defect_dir;

	for(a=0;a<num_swaps;a++)
		{
		b=a+1;
		metro_swap_prob = 0.0;
		// compute action difference between replica a and b
		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) reduction(+:metro_swap_prob) private(n,aux_i,j)
		#endif
		for(n=0;n<(swap_rectangle->d_vol_rect);n++)
			{
			long r = swap_rectangle->rect_sites[n]; // action changes only in the first neighborhood of the defect
			for(aux_i=0; aux_i<STDIM-1; aux_i++)
				{
				j = perp_dir[param->d_defect_dir][aux_i];
				metro_swap_prob += delta_action_swap(GC, geo, param, r, i, j, a, b);
				}
			}
	
		// do the swap
		metro_swap_prob=exp(-metro_swap_prob); // metropolis swap probability
		metropolis_single_swap(GC,a,b,metro_swap_prob,acc_counters);
		}
	}
*/

// metropolis step to swap replica a and b with probability p, including the twist factors
void metropolis_single_swap(Gauge_Conf *GC, int const a, int const b, double const p, Acc_Utils *acc_counters)
	{
	// acceptance initialized to 1
	int acc=1;
	// increase counter of tried swaps for replicas (a, a+1)
	acc_counters->num_swap[a]++;

	// Metropolis test: if p<1 => acc=1 with probability p, if p>=1 acc=1 (already assigned)
	if(p<1)
	{
		double random_number=casuale();
		if(random_number>p)
		{
			acc=0;
		}
	}
	
	// if Metropolis is accepted, swap replicas, including the twist factors
	if(acc==1)
		{
		// swap of configurations
		GAUGE_GROUP **aux;
		double complex **aux_Z;
		aux=GC[a].lattice;
		aux_Z=GC[a].Z;
		GC[a].lattice=GC[b].lattice;
		GC[a].Z=GC[b].Z;
		GC[b].lattice=aux;
		GC[b].Z=aux_Z;
		acc_counters->num_accepted_swap[a]++; // increase counter of successfull swaps for replicas (a, a+1)

		// swap of labels
		int aux_label;
		aux_label=GC[a].conf_label;
		GC[a].conf_label=GC[b].conf_label;
		GC[b].conf_label=aux_label;
		}
	}

// translation of one lattice spacing of the configuration, including the twist factors
// direction is chosen randomly, verse is always positive
void conf_translation(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param)
	{
	double aux;
	int dir,i;
	long s;
	Gauge_Conf aux_conf;
  
	// extract random direction
	aux=STDIM*casuale();
	for(i=0;i<STDIM;i++)
	{
		if ( (aux>=i) && (aux<(i+1)) ) dir=i;
	}

	// copy the conf in an auxiliary one (should be defined outside and passed to the function?), including the twist factors
	init_gauge_conf_from_gauge_conf(&aux_conf, GC, param); // now aux_conf=GC

	// translation in direction +dir, including the twist factors
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s)
	#endif
	for(s=0; s<(param->d_n_planes)*(param->d_volume); s++)
	{
		// s = j * volume + r
		long r = s % (param->d_volume);
		int j = (int) ( (s-r)/(param->d_volume) );
		if(j<STDIM) equal(&(GC->lattice[r][j]), &(aux_conf.lattice[nnm(geo,r,dir)][j]) );
		GC->Z[r][j] = aux_conf.Z[nnm(geo,r,dir)][j];
	}

	// free auxiliary conf, including the twist factors
	free_gauge_conf(&aux_conf, param);
	free_twist_cond(&aux_conf, param);
	}
	
void init_swap_acc_arrays(Acc_Utils *acc_counters, GParam const * const param)
  {
	if(param->d_N_replica_pt==1)
		{
		acc_counters->num_accepted_swap=NULL;
		acc_counters->num_swap=NULL;
		}
	else
		{
		int i,err;
	
		err=posix_memalign( (void **) &(acc_counters->num_accepted_swap), (size_t) INT_ALIGN, (size_t) ((param->d_N_replica_pt)-1) * sizeof(long));
		if(err!=0)
			{
			fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}

		err=posix_memalign( (void **) &(acc_counters->num_swap), (size_t) INT_ALIGN, (size_t) ((param->d_N_replica_pt)-1) * sizeof(long));
		if(err!=0)
			{
			fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	
		for(i=0;i<((param->d_N_replica_pt)-1);i++) 
			{
			acc_counters->num_accepted_swap[i]=0;
			acc_counters->num_swap[i]=0;
			}
		}
  }
	
void end_swap_acc_arrays(Acc_Utils *acc_counters, GParam const * const param)
	{
	if(param->d_N_replica_pt>1)
		{
		free(acc_counters->num_accepted_swap);
		free(acc_counters->num_swap);
		}
	else
		{
		(void) acc_counters; // to suppress compiler warning of unused variable
		}
	}
  
void print_acceptances(Acc_Utils const * const acc_counters, GParam const * const param)
  {
	if(param->d_N_replica_pt==1)
		{
		(void) acc_counters; // to suppress compiler warning of unused variable
		}
	else
		{
		FILE *fp;	
		double acc,err_acc;
		int r;
  
		fp=fopen(param->d_swap_acc_file, "w");
		fprintf(fp, "#swap_from    swap_to    c_1    c_2    acceptance(%%)    err_acceptance(%%)    swap_accepted    swap_tried\n");
		for(r=0;r<((param->d_N_replica_pt)-1);r++)
			{
			if(acc_counters->num_swap[r]!=0)
				{
				acc = ( (double) (acc_counters->num_accepted_swap[r]) ) / ( (double) (acc_counters->num_swap[r]) ) ;
				err_acc = sqrt( acc * (1.0-acc) / ( ( (double) (acc_counters->num_swap[r]) ) - 1.0 ) );
				}
			else
				{
				acc=0.0;
				err_acc=0.0;
				}
			fprintf(fp,"%d    %d    %lf    %lf    %lf    %lf    %ld    %ld\n", r, r+1, param->d_defect_beta[r], param->d_defect_beta[r+1], acc*100.0, err_acc*100.0, acc_counters->num_accepted_swap[r], acc_counters->num_swap[r]);
			}
		fclose(fp);	  
		}
	}

void init_swap_track_file(FILE **swaptrackfilep, GParam const * const param)
{
	if (param->d_N_replica_pt > 1)
	{
		if (param->d_start==2) // starting run from saved conf
		{
			*swaptrackfilep=fopen(param->d_swap_tracking_file, "r");
			if(*swaptrackfilep!=NULL) // file exists -> close it and re-open it in append mode
			{
				fclose(*swaptrackfilep);
				*swaptrackfilep=fopen(param->d_swap_tracking_file, "a");
			}
			else // file does not exist -> create it and write first line
			{
				*swaptrackfilep=fopen(param->d_swap_tracking_file, "w");
				fprintf(*swaptrackfilep, "# MC_step    conf_labels\n");
				fflush(*swaptrackfilep);
			}
		}
		else // starting run from scratch
		{
			*swaptrackfilep=fopen(param->d_swap_tracking_file, "w");
			fprintf(*swaptrackfilep, "# MC_step    conf_labels\n");
			fflush(*swaptrackfilep);
		}
	}
	else // no need of this file if num_replica = 1
	{
		(void) swaptrackfilep; // to suppress compiler warning of unused variable
		(void) param; // to suppress compiler warning of unused variable
	}
}
  
void print_conf_labels(FILE *fp, Gauge_Conf const * const GC, GParam const * const param)
{
	if (param->d_N_replica_pt>1)
	{
		fprintf(fp, "%ld      ",GC[0].update_index);
		for(int r=0;r<(param->d_N_replica_pt);r++) fprintf(fp,"%d ",GC[r].conf_label);
		fprintf(fp,"\n");
		fflush(fp);
	}
	else
	{
		(void) fp; // to suppress compiler warning of unused variable
		(void) GC; // to suppress compiler warning of unused variable
		(void) param;  // to suppress compiler warning of unused variable
	}
}
#endif
