#ifndef GAUGE_CONF_UPDATE_MULTICANONIC_C
#define GAUGE_CONF_UPDATE_MULTICANONIC_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// read topo potential from file
void read_topo_potential(double **grid, GParam const * const param)
{
	int i, j, err;
	double x, V;
	FILE *fp;

	fp=fopen(param->d_topo_potential_file, "r");
	if( fp==NULL )
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_topo_potential_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}	
	
	err=posix_memalign((void **) grid, (size_t) DOUBLE_ALIGN, (size_t) param->d_n_grid * sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the grid array!\n");
		exit(EXIT_FAILURE);
	}

	// read x and V(x) from topo_potential file
	for (i=0; i<param->d_n_grid; i++)
	{
		err=fscanf(fp, "%lf %lf", &x, &V); 
		if(err!=2)
		{
			printf("Error: can't read the %d-th element of the file (%s, %d)\n", i,__FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		j=(int)floor((x+param->d_grid_max+(param->d_grid_step/2.0))/param->d_grid_step);
		if (i!= j)
		{
			printf("Error: found %d (%lf) when expecting %d (%s, %d)\n", j, x, i, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		(*grid)[i]=V;
	}
	fclose(fp);
}

// compute topo potential in point x
double compute_topo_potential(double const x, double const * const grid, GParam const * const param)
{	  
	int i_grid;
	double x0, m, q;

	// find index of nearest grid point to x
	i_grid=(int)(floor((x+param->d_grid_max)/param->d_grid_step));

	if(i_grid>=0 && i_grid<param->d_n_grid) // if x inside the barriers compute V(x) with a linear interpolation
	{
		// perform linear interpolation: V(x) = V(x0) + [dV/dx|(x0)] (x-x0)
		x0=i_grid*param->d_grid_step-param->d_grid_max;
		m=(grid[i_grid+1]-grid[i_grid])/param->d_grid_step; // dV/dx|(x0) = [ V(x0+step) - V(x0) ] / step
		q=grid[i_grid]-m*x0; // V(x0) - [dV/dx|(x0)] x0
		return q+m*x; // V(x0) + [dV/dx|(x0)] (x-x0) 
	}
	else // if x outside the barriers just saturate to extreme values
	{
		if(i_grid<0) return grid[0];
		else return grid[param->d_n_grid-1];
	}
}

// initialize multicanonic acceptances file
void init_multicanonic_acc_file(FILE **filep, GParam const * const param)
{
	int i;
	
	(*filep)=fopen(param->d_multicanonic_acc_file, "r");
	if((*filep)!=NULL) // file exists
	{
		fclose((*filep));
		(*filep)=fopen(param->d_multicanonic_acc_file, "a");
	}
	else // file doesn't exist, write first line
	{
		(*filep)=fopen(param->d_multicanonic_acc_file, "w");
		fprintf((*filep), "# %f ", param->d_beta);
		for(i=0; i<STDIM; i++) fprintf(*filep, "%d ", param->d_size[i]);
		fprintf((*filep), "\n");
	}
	fflush(*filep);
}

void init_multicanonic_acc_arrays(Acc_Utils *acc_counters, GParam const * const param)
{
		int i,err;
	
		err=posix_memalign( (void **) &(acc_counters->num_accepted_metro_multicanonic), (size_t) INT_ALIGN, (size_t) (param->d_N_replica_pt) * sizeof(long));
		if(err!=0)
			{
			fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}

		err=posix_memalign( (void **) &(acc_counters->num_metro_multicanonic), (size_t) INT_ALIGN, (size_t) (param->d_N_replica_pt) * sizeof(long));
		if(err!=0)
			{
			fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	
		for(i=0;i<(param->d_N_replica_pt);i++) 
			{
			acc_counters->num_accepted_metro_multicanonic[i]=0;
			acc_counters->num_metro_multicanonic[i]=0;
			}
}
			
void end_multicanonic_acc_arrays(Acc_Utils *acc_counters)
{
	free(acc_counters->num_accepted_metro_multicanonic);
	free(acc_counters->num_metro_multicanonic);
}

// print metropolis acceptance of multicanonic algorithm on file
void print_multicanonic_acceptance(Gauge_Conf const * const GC, GParam const * const param, Acc_Utils const * const acc_counters, FILE *multicanonic_acc_filep)
{
	double mean_acc;
	if(GC->update_index % param->d_measevery == 0)
	{
		fprintf(multicanonic_acc_filep, "%ld ", GC->update_index);
		for(int i=0; i<param->d_N_replica_pt; i++)
		{
			if (acc_counters->num_metro_multicanonic[i] != 0) mean_acc = ( (double) acc_counters->num_accepted_metro_multicanonic[i] ) / ( (double) acc_counters->num_metro_multicanonic[i] );
			else mean_acc = 0.0;
			fprintf(multicanonic_acc_filep, "%lf ", 100.0*mean_acc);
		}
		fprintf(multicanonic_acc_filep, "\n");
		fflush(multicanonic_acc_filep);
	}
}

// perform a single step of parallel tempering with hierarchic update with a multicanonical approach
void multicanonic_parallel_tempering_with_hierarchical_update(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
												 Rectangle const * const most_update, Rectangle const * const clover_rectangle,
												 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters, double const * const grid, FILE *acc_filep)												 
{
	// first hierarc level is 0
	int start_hierarc=0;
	
	// set multicanonic Metropolis acceptance counters to zero to compute mean acc over single updating step
	for(int i=0;i<param->d_N_replica_pt; i++)
	{
		acc_counters->num_accepted_metro_multicanonic[i] = 0;
		acc_counters->num_metro_multicanonic[i] = 0;
	}
	
	// Parallel tempering updating step: full update + hierarchical update + swaps and translations after every sweep for every replica
	multicanonic_update_with_defect(GC, geo, param, grid, acc_counters); // full update of all replicas
	if(param->d_N_replica_pt>1)
	{
		swap(GC, geo, param, swap_rectangle, acc_counters); // swap all replicas
		conf_translation(&(GC[0]), geo, param); // translation of periodic replica (GC[0])
		if(param->d_N_hierarc_levels>0)
			multicanonic_hierarchical_update_rectangle_with_defect(GC, geo, param, start_hierarc, most_update, clover_rectangle,
																															swap_rectangle, acc_counters, grid); // hierarchic update
	}
	
	// increase update index of all replicas
	for(int i=0;i<param->d_N_replica_pt; i++) GC[i].update_index++;
	
	// print mean multicanonic acceptance over a single updating step
	print_multicanonic_acceptance(GC, param, acc_counters, acc_filep);
}

// update all replica only on a given rectangle in the presence of a defect (MODIFICARE)
void multicanonic_update_with_defect(Gauge_Conf * GC, Geometry const * const geo, GParam const * const param, double const * const grid, Acc_Utils *acc_counters)
{
for(int i=0; i<STDIM; i++)
      {
      if(param->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long s;
   int j, dir, err;
	 int num_replica = param->d_N_replica_pt; // just an auxiliary variable
	 long *sum_acc, *count_metro;

	 // init aux variables to compute mean multicanonic acc
	 	err=posix_memalign( (void **) &(sum_acc), (size_t) INT_ALIGN, (size_t) num_replica * sizeof(long));
	if(err!=0)
		{
		fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
			err=posix_memalign( (void **) &(count_metro), (size_t) INT_ALIGN, (size_t) num_replica * sizeof(long));
	if(err!=0)
		{
		fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
		
		for(int i=0; i<param->d_N_replica_pt; i++)
		{
			sum_acc[i]=0;
			count_metro[i]=0;
		}

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
				compute_clovers_replica(GC, geo, param, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
      #endif
			for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
				{
				// s = i * volume/2 + r
				long r = s % ( (param->d_volume)/2 ); // site index
				int i = (int) ( (s-r) / ( (param->d_volume)/2 ) ); // replica index
				int acc_metro;
				acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir, grid);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
      #endif 
			for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
				{
				// s = i * volume/2 + aux ; aux = r - volume/2
				long aux = s % ( (param->d_volume)/2 );
				long r = (param->d_volume/2) + aux; // site index
				int i = (int) ( (s-aux) / ( (param->d_volume)/2 ) ); // replica index
				int acc_metro;
				acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir, grid);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
				compute_clovers_replica(GC, geo, param, dir);
      #endif

      for(j=0; j<param->d_overrelax; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
         #endif 
				 for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
						{
						// s = i * volume/2 + r
						long r = s % ( (param->d_volume)/2 ); // site index
						int i = (int) ( (s-r) / ( (param->d_volume)/2 ) ); // replica index
						int acc_metro;
            acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir, grid);
						sum_acc[i] += acc_metro;
						count_metro[i]++;
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
         #endif 
				 for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)/2); s++)
						{
						// s = i * volume/2 + aux ; aux = r - volume/2
						long aux = s % ( (param->d_volume)/2 );
						long r = (param->d_volume/2) + aux; // site index
						int i = (int) ( (s-aux) / ( (param->d_volume)/2 ) ); // replica index
						int acc_metro;
            acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir, grid);
						sum_acc[i] += acc_metro;
						count_metro[i]++;
            }
         }
      }
			
		// add number of accepted and number of proposed Metropolis multicanonic tests
		for(int i=0; i<param->d_N_replica_pt; i++)
		{
				acc_counters->num_accepted_metro_multicanonic[i] += sum_acc[i];
        acc_counters->num_metro_multicanonic[i] += count_metro[i];
		}
		
		// free aux arrays
		free(sum_acc);
		free(count_metro);
   
   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(s, dir)
   #endif 
   for(s=0; s<((param->d_N_replica_pt)*(param->d_volume)); s++)
      {
			// s = i * volume + r
			long r = s % param->d_volume;
			int i = (int) ( (s-r) / (param->d_volume) );
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC[i].lattice[r][dir]));
         } 
      }

   }

// hierarchical update functions

// update all replica only on a given rectangle in the presence of a defect (MODIFICARE)
void multicanonic_update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
																								Rectangle const * const most_update, Rectangle const * const clover_rectangle, double const * const grid,Acc_Utils *acc_counters)
	{
	
	long s, num_even, num_odd;
	int j,dir,err;
	 int num_replica = param->d_N_replica_pt; // just an auxiliary variable
	 long *sum_acc, *count_metro;

	 // init aux variables to compute mean multicanonic acc
	 	err=posix_memalign( (void **) &(sum_acc), (size_t) INT_ALIGN, (size_t) num_replica * sizeof(long));
	if(err!=0)
		{
		fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
			err=posix_memalign( (void **) &(count_metro), (size_t) INT_ALIGN, (size_t) num_replica * sizeof(long));
	if(err!=0)
		{
		fprintf(stderr, "Problems in allocating the acceptances array (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
		
		for(int i=0; i<param->d_N_replica_pt; i++)
		{
			sum_acc[i]=0;
			count_metro[i]=0;
		}
	
	#ifndef THETA_MODE
		(void) clover_rectangle; // to avoid compiler warning of unused variable
	#endif	
		
	/* Check if there's at least one even dimension of the rectangle, i.e. check if d_vol_rect is even.
	   If there's at least one even dimension: d_vol_rect/2 even sites and d_vol_rect/2 odd sites.
		 Otherwise: (d_vol_rect+1)/2 even sites and (d_vol_rect-1)/2 odd sites. */

	long is_even = ( most_update->d_vol_rect ) % 2;

	num_even = ( most_update->d_vol_rect + is_even ) / 2; // number of even sites
	num_odd  = ( most_update->d_vol_rect - is_even ) / 2; // number of odd sites

	// heatbath
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
			compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif 
		for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_even + n
			long n = s % num_even;               // site index on rectangle
			long r = most_update->rect_sites[n]; // site index on lattice
			int i = (int) ( (s-n) / num_even );  // replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir, grid);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			}

		#ifdef OPENMP_MODE
		#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
		#endif 
		for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
			{
			// s = i * num_odd + aux; aux = n - num_even
			long aux = s % num_odd;
			long n = aux + num_even;             // site index on rectangle
			long r = most_update->rect_sites[n]; // site index on lattice
			int i = (int) ( (s-aux) / num_odd ); // replica index
			int acc_metro;
			acc_metro = multicanonic_heatbath_with_defect(&(GC[i]), geo, param, r, dir, grid);
			sum_acc[i] += acc_metro;
			count_metro[i]++;
			} 
		}

	// overrelax
	for(dir=0; dir<STDIM; dir++)
		{
		#ifdef THETA_MODE
			compute_clovers_replica_rect(GC, geo, param, dir, clover_rectangle);
		#endif

		for(j=0; j<param->d_overrelax; j++)
			{
			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<(num_even*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_even + n
				long n = s % num_even;               // site index on rectangle
				long r = most_update->rect_sites[n]; // site index on lattice
				int i = (int) ( (s-n) / num_even );  // replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir, grid);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}

			#ifdef OPENMP_MODE
			#pragma omp parallel for num_threads(NTHREADS) private(s) reduction(+:sum_acc[:num_replica]) reduction(+:count_metro[:num_replica])
			#endif 
			for(s=0; s<(num_odd*(param->d_N_replica_pt)); s++)
				{
				// s = i * num_odd + aux; aux = n - num_even
				long aux = s % num_odd; 
				long n = aux + num_even;             // site index on rectangle
				long r = most_update->rect_sites[n]; // site index on lattice
				int i = (int) ( (s-aux) / num_odd ); // replica index
				int acc_metro;
				acc_metro = multicanonic_overrelaxation_with_defect(&(GC[i]), geo, param, r, dir, grid);
				sum_acc[i] += acc_metro;
				count_metro[i]++;
				}
			}
		}
 
		// add number of accepted and number of proposed Metropolis multicanonic tests
		for(int i=0; i<param->d_N_replica_pt; i++)
		{
				acc_counters->num_accepted_metro_multicanonic[i] += sum_acc[i];
        acc_counters->num_metro_multicanonic[i] += count_metro[i];
		}
		
		// free aux arrays
		free(sum_acc);
		free(count_metro);
 
	// final unitarization
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s, dir)
	#endif 
	for(s=0; s<((most_update->d_vol_rect)*(param->d_N_replica_pt)); s++)
		{
			for(dir=0; dir<STDIM; dir++)
				{
				// s = i * volume_rect + n
				long n = s % (most_update->d_vol_rect); // site index on rectangle
				long r = most_update->rect_sites[n];    // site index on lattice
				int i = (int) ( (s-n) / (most_update->d_vol_rect) ); // replica index
				unitarize(&(GC[i].lattice[r][dir]));
				} 
		}

	}

// perform a hierarchical update on all rectangles
void multicanonic_hierarchical_update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
														 int const hierarc_level, 
														 Rectangle const * const most_update,
														 Rectangle const * const clover_rectangle,
														 Rectangle const * const swap_rectangle,
														 Acc_Utils *acc_counters,
														 double const * const grid)
	{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
			{
	    multicanonic_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]),
																								grid, acc_counters);
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
	    conf_translation(&(GC[0]), geo, param);
	    }
		} // end if
	else
		{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
			{
			multicanonic_update_rectangle_with_defect(GC,geo,param, &(most_update[hierarc_level]), &(clover_rectangle[hierarc_level]),
														grid, acc_counters);   
			if(param->d_N_replica_pt>1) swap(GC, geo, param, swap_rectangle, acc_counters);
			conf_translation(&(GC[0]), geo, param);
			multicanonic_hierarchical_update_rectangle_with_defect(GC,geo,param,(hierarc_level+1),most_update,clover_rectangle,
																	swap_rectangle,acc_counters, grid);
			}
		} // end else
	}

// perform an update with overrelaxation and multicanonic Metropolis step
int multicanonic_overrelaxation_with_defect(Gauge_Conf *GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i, double const * const grid)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  (void) param; // just to avoid warnings

  GAUGE_GROUP stap;

  #ifndef THETA_MODE
    calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
  #else
    calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
  #endif
  
  // store old link before update
  GAUGE_GROUP old_link = GC->lattice[r][i];
	
  // perform usual single-link over-relaxation update
  single_overrelaxation(&(GC->lattice[r][i]), &stap);
  
  // accept/reject new link with multicanonic Metropolis step
  int acc = multicanonic_Metropolis_step_single_link(GC, geo, param, r, i, grid, old_link);
	return acc;
  }
 
// perform an update with heatbath in the presence of a defect with multicanonic Metropolis step
int multicanonic_heatbath_with_defect(Gauge_Conf *GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i, double const * const grid)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  GAUGE_GROUP stap;

  #ifndef THETA_MODE
    calcstaples_wilson_with_defect(GC, geo, param, r, i, &stap);
  #else
    calcstaples_with_topo_with_defect(GC, geo, param, r, i, &stap);
  #endif

  // store old link before update
  GAUGE_GROUP old_link = GC->lattice[r][i];

  // perform usual single-link over-heat-bath update
  single_heatbath(&(GC->lattice[r][i]), &stap, param);
  
  // accept/reject new link with multicanonic Metropolis step
  int acc = multicanonic_Metropolis_step_single_link(GC, geo, param, r, i, grid, old_link);
	return acc;
  }

// initialize topo_charge for all replicas
void init_topo_charge(Gauge_Conf * GC, Geometry const * const geo, GParam const * const param)
{
	for (int i=0; i<param->d_N_replica_pt; i++)
		refresh_topo_charge(&(GC[i]), geo, param);
}

// refresh topo_charge for non-periodic replica
void refresh_topo_charge_replica(Gauge_Conf * GC, Geometry const * const geo, GParam const * const param)
{
	// starts from 1 because the topo charge of replica with i=0 is refreshed before measuring observables
	for (int i=1; i<param->d_N_replica_pt; i++)
		refresh_topo_charge(&(GC[i]), geo, param);
}

// refresh topo_charge of given replica
void refresh_topo_charge(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param)
{
	double Q;
	Q = topcharge(GC, geo, param);
	GC->stored_topo_charge = Q;
}

// compute single clover insertion and store it in <clover_insertion>
void compute_single_clover_insertion(GAUGE_GROUP * clover_insertion, Gauge_Conf const * const GC,
                     Geometry const * const geo, GParam const * const param, long r, int i, int j)
{
	GAUGE_GROUP aux;
	clover(GC, geo, param, r, i, j, &aux); // aux = clover[r][i][j]
	equal(clover_insertion, &aux); // clover_insertion = aux
	minus_equal_dag(clover_insertion, &aux);  // clover_insertion -= aux^{dag}
	// ==> clover_insertion = clover[r][i][j] - clover[r][i][j]^dag = 2i Im{ clover[r][i][j] }
}

// compute staple of topological charge relative to link (r,i)
void compute_topostaple_alone(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, long r, int i, GAUGE_GROUP * topo_stap)
{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if(i >= STDIM)
	{
		fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	#endif

	if(STDIM!=4)
	{
		fprintf(stderr, "Error: topological charge can be used only in 4 dimensions! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	GAUGE_GROUP link1, link2, link3, link12, stap, aux, clover_insertion;
	const double coeff = 1.0/(128.0*PI*PI);
	long k;
	int j, l;
	int i0, j0;
	int sood1[4][4], sood2[4][4]; // Signed Ordered Orthogonal Directions (SOOD)

	zero(topo_stap); // topo_stap=0

	// the clover topological charge is written as
	// -1/(128 pi^2) \sum_{ind. perm.} ReTr(Q_{\mu\nu}(Q-Q^{dag})_{sood1[\mu][\nu] sood2[\mu][\nu]} )
	// where Q_{\mu\nu} here stands for the clover on plane (\mu \nu)
	// the independent permutations are 3, here we use: 0123 0231 0312

	sood1[0][1] = 2;
	sood2[0][1] = 3;
	sood1[1][0] = 3;
	sood2[1][0] = 2;

	sood1[0][2] = 3;
	sood2[0][2] = 1;
	sood1[2][0] = 1;
	sood2[2][0] = 3;

	sood1[0][3] = 1;
	sood2[0][3] = 2;
	sood1[3][0] = 2;
	sood2[3][0] = 1;

	sood1[1][2] = 0;
	sood2[1][2] = 3;
	sood1[2][1] = 3;
	sood2[2][1] = 0;

	sood1[1][3] = 2;
	sood2[1][3] = 0;
	sood1[3][1] = 0;
	sood2[3][1] = 2;

	sood1[2][3] = 0;
	sood2[2][3] = 1;
	sood1[3][2] = 1;
	sood2[3][2] = 0;

	for(l=i+1; l< i + STDIM; l++)
	{
		j = (l % STDIM);

		i0=sood1[i][j];
		j0=sood2[i][j];

//
//       i ^
//         |   (1)
//     (b) +----->-----+ (c)
//         |           |
//         |           |
//         |           V (2)
//         |           |
//         |           |
//     (a) +-----<-----+-->   j
//       r     (3)    (d)
//
			
		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
		equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
		equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

		times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
		times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

		// clover insertion in (a)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, r, i0, j0);
		times(&aux, &stap, &clover_insertion); // stap*clover
		plus_equal(topo_stap, &aux);

		// clover insertion in (b)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, i), i0, j0);
		times(&aux, &clover_insertion, &stap);  // clover*stap
		plus_equal(topo_stap, &aux);

		// clover insertion in (c)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, nnp(geo, r, i), j), i0, j0);
		times(&aux, &link1, &clover_insertion);  // link1*clover
		times_equal_dag(&aux, &link2);       // *=link2^{dag}
		times_equal_dag(&aux, &link3);       // *=link3^{dag}
		plus_equal(topo_stap, &aux);

		// clover insertion in (d)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, j), i0, j0);
		times(&aux, &link12, &clover_insertion);  // link1*link2*quadri
		times_equal_dag(&aux, &link3);          // *=link3^{dag}
		plus_equal(topo_stap, &aux);

//
//       i ^
//         |   (1)
//     (d) +----<------+ (a)
//         |           |
//         |           |
//     (2) V           |
//         |           |
//         |           | (b)
//     (c) +------>----+--->j
//        k     (3)    r
//

		k=nnm(geo, r, j); // k = r - j

		// non-topo staple
		equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
		equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
		equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

		times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
		times(&stap, &link12, &link3);        // stap=link12*link3

		// clover insertion in (a)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, r, i), i0, j0);
		times(&aux, &clover_insertion, &stap); // clover*stap
		minus_equal(topo_stap, &aux);

		// clover insertion in (b)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, r, i0, j0);
		times(&aux, &stap, &clover_insertion); // stap*clover
		minus_equal(topo_stap, &aux);

		// clover insertion in (c)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, k, i0, j0);
		times(&aux, &link12, &clover_insertion); // link1^{dag}*link2^{dag}*clover
		times_equal(&aux, &link3);                            // *=link3
		minus_equal(topo_stap, &aux);

		// clover insertion in (d)
		compute_single_clover_insertion(&clover_insertion, GC, geo, param, nnp(geo, k, i), i0, j0);
		times_dag1(&aux, &link1, &clover_insertion);  // link1^{dag}*clover
		times_equal_dag(&aux, &link2);            // *=link2^{dag}
		times_equal(&aux, &link3);                // *=link3
		minus_equal(topo_stap, &aux);
	}
	times_equal_real(topo_stap, coeff); // topo_stap *= coeff
}

// compute the variation of the clover topological charge when the link (r,i) is updated starting from <old_link>
double delta_Q_upd(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, long r, int i, GAUGE_GROUP old_link)
{
	// compute delta_Q
	GAUGE_GROUP q, topo_stap;
	double delta_Q = 0.0;
	
	// compute topo staple and store it in topo_stap
	compute_topostaple_alone(GC, geo, param, r, i, &topo_stap);
	
	// compute contribution of new link
	times(&q, &topo_stap, &(GC->lattice[r][i])); // q = (topo_stap * new_link) 
	delta_Q += retr(&q)*((double) NCOLOR); // delta_Q += [ retr(topo_stap * new_link) / NCOLOR ] * NCOLOR (retr automatically adds a 1/NCOLOR factor)
	
	// compute contribution of old_link
	times(&q, &topo_stap, &old_link); // q = (topo_stap * old_link)
	delta_Q -= retr(&q)*((double) NCOLOR); // delta_Q -= [ retr(topo_stap * old_link) / NCOLOR ] * NCOLOR
	
	return delta_Q;
}

// compute the multicanonic Metropolis probability p=exp(delta V) where V=topo potential
double metropolis_prob_multicanonic(double const Q_new, double const Q_old, GParam const * const param, double const * const grid)
{
	double V_old, V_new, delta_V, p;

	V_old = compute_topo_potential(Q_old, grid, param);
	V_new = compute_topo_potential(Q_new, grid, param);
	delta_V = V_new-V_old;
	p = exp(-delta_V);
	return p;
}

// TO BE MODIFIED
//  Metropolis test with p_metro=exp(- delta topo_potential)
int multicanonic_Metropolis_step_single_link(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
												long r, int i, double const * const grid, GAUGE_GROUP old_link)
{
	// perform multicanonic Metropolis test
	double Q_old = GC->stored_topo_charge;
	double Q_new = Q_old + delta_Q_upd(GC, geo, param, r, i, old_link);
	double p = metropolis_prob_multicanonic(Q_new, Q_old, param, grid);
	
	// Metropolis test: p < 1 --> acc=1 with probability p, else --> acc=1
	int acc = 1;
	if(p < 1)
	{
		double random_number=casuale();
		if(random_number > p) acc = 0;
	}

	if (acc == 0) GC->lattice[r][i] = old_link;	// if Metropolis is refused go back to original link
	else GC->stored_topo_charge = Q_new;	      // if Metropolis is accepted store the new topological charge of the conf

	return acc;
}

#endif
