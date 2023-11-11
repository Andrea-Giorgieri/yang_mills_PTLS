#ifndef DEBUG_TWIST_C
#define DEBUG_TWIST_C

#include"../config.h"
#include"../include/macro.h"
#include"../include/gparam.h"
#include"../include/random.h"
#include"../include/function_pointers.h"
#include"../include/sun.h"
#include"../include/gauge_conf.h"

#include<math.h>
#include<stdlib.h>

void printfc(double complex z);
void conf_translation_dir(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param, int dir); //dir not random, debug only

int main(void)
	{
	double plaqs, plaqt, plaqs_new, plaqt_new, C_a, C_b, energy, action, swapped_action, re_tr_plaq_a, re_tr_plaq_b, K_a, K_b;
	double complex z, w, trace_calcstaples, trace_plaquettep, trace_plaquettep_swap, trace_clover, Z_a, Z_b;
	char in_file[] = "input_file_ym_pt";
	int i, j, k, err=0;
	long r=0;
	
	Gauge_Conf *GC;
	GAUGE_GROUP M, N;
	Geometry geo;
	GParam param;
	Rectangle swap_rectangle;
	Rectangle *most_update, *clover_rectangle;
	Acc_Utils acc_counters;
	int L_R_swap=1;
	
	readinput(in_file, &param);
	if(param.d_N_replica_pt<2) param.d_N_replica_pt = 2;
	param.d_plaquette_meas=1;
	param.d_clover_energy_meas=1;
	param.d_charge_meas=1;
	param.d_chi_prime_meas = 1;
	param.d_polyakov_meas = 1;
	param.d_topcharge_tcorr_meas = 1;
	
	initrand(param.d_randseed);
	init_indexing_lexeo();
	init_geometry(&geo, &param);
	init_gauge_conf_replica(&GC, &param);
	init_rect_hierarc(&most_update, &clover_rectangle, &param);
	init_rect(&swap_rectangle, L_R_swap, &param);
	init_swap_acc_arrays(&acc_counters, &param);
		
	printf("\n*******************************\n");
	printf("PROGRAM FOR THE DEBUG OF TWIST\n");
	printf("*******************************\n\n");

	printf("N=%s", QUOTEME(NCOLOR));
	
	printf("\n\n");
	printf("PRINT TWIST FACTORS AT ORIGIN ...\n\n");
	
	for(i=0; i<STDIM; i++)
	{
		for(j=i+1; j<STDIM; j++) 
		{
			z = GC[0].Z[r][dirs_to_si(i,j)];
			printf("Z(%d, %d) = %g%s%gi\n", i, j, creal(z), (cimag(z)>=0.0f)? "+":"", cimag(z));
			z = GC[0].Z[r][dirs_to_si(j,i)];
			printf("Z(%d, %d) = %g%s%gi\n", j, i, creal(z), (cimag(z)>=0.0f)? "+":"", cimag(z));
		}
	}
	
	printf("\n\n");
	printf("VERIFY THAT TWISTED COLD START GIVES ZERO ENERGY ...\n\n");
	
	param.d_start = 3;
	init_gauge_conf_replica(&GC, &param);
	for(r=0; r<(param.d_volume); r++)
	{
		energy = 0;
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++)
			{
				clover(&(GC[0]), &geo, &param, r, i, j, &M);
				ta(&M);
				equal(&N, &M);
				times_equal(&M, &N);
				energy += -NCOLOR*retr(&M)/16.0;
			}
		if (fabs(energy) > MIN_VALUE)
		{
			err = 1;
			printf("	ERROR!!!!!!!!!!!\n\n");
			printf("	Pos %ld Energy = %g\n\n", r, energy);
		}
	}
	if (err == 0) printf("Vacuum Energy :		OK\n");
	err = 0;
	r=0;
	
	for(r=0; r<(param.d_volume); r++)
	{
		energy = 0;
		for(i=0; i<STDIM; i++)
			for(j=i+1; j<STDIM; j++) 
			{
				trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r, i, j);
				if (cabs(trace_plaquettep - 1.0) > MIN_VALUE)
				{
					err = 1;
					printf("	ERROR!!!!!!!!!!!\n\n");
					printf("	Pos %ld dirs %d %d Plaquette = ", r, i, j);
					printfc(trace_plaquettep);
					printf("\n\n");
				}
			}
	}
	if (err == 0) printf("Vacuum Plaquette :	OK\n");
	err = 0;
	r=0;
	
	printf("\n\n");	
	printf("VERIFY THAT plaquettep GETS CONJUGATED SWAPPING THE DIRS OF A PLAQUETTE AT ORIGIN...\n\n");
	
	param.d_start = 1;
	init_gauge_conf_replica(&GC, &param);
	for(i=0; i<STDIM; i++)
	{
		for(j=i+1; j<STDIM; j++)
		{
			trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r, i, j);
			trace_plaquettep_swap = plaquettep_complex(&(GC[0]), &geo, &param, r, j, i);
			
			if(cabs(trace_plaquettep_swap-conj(trace_plaquettep))<MIN_VALUE)
			{
				printf("Trace plane (%d, %d) :	OK\n", i, j);
			}
			else
			{
				printf("	ERROR!!!!!!!!!!!\n\n");
				printf("	plane (%d, %d) plaquettep = ", i, j);
				printfc(trace_plaquettep);
				printf(", plaquettep_swap = ");
				printfc(trace_plaquettep_swap);
				printf("\n\n");
			}
		}
	}
	
	printf("\n\n");
	printf("VERIFY THAT calcstaples AND plaquettep GIVE THE SAME RESULT FOR THE TRACE OF A LINK AT ORIGIN ...\n\n");
	
	for(i=0; i<STDIM; i++)
	{
		calcstaples_wilson(&(GC[0]), &geo, &param, r, i, &M);
		times_equal(&M, &(GC[0].lattice[r][i]));
		trace_calcstaples = retr(&M) + I*imtr(&M);
		
		trace_plaquettep = 0.0 + I*0.0;
		for(j=0; j<STDIM; j++)
		{
			if(j != i)
			{
				trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, r, j, i);
				trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r, j), i, j);
			}
		}
		
		if(cabs(trace_calcstaples-trace_plaquettep)<MIN_VALUE)
		{
			printf("Trace dir %d :	OK\n", i);
		}
		else
		{
			printf("	ERROR!!!!!!!!!!!\n\n");
			printf("	dir %d calcstaples = ", i);
			printfc(trace_calcstaples);
			printf(", plaquettep = ");
			printfc(trace_plaquettep);
			printf("\n\n");
		}
	}
	
	printf("\n\n");
	printf("VERIFY THAT clover AND plaquettep GIVE THE SAME RESULT FOR THE TRACE OF A CLOVER AT ORIGIN ...\n\n");
	
	for(i=0; i<STDIM; i++)
	{
		for(j=i+1; j<STDIM; j++)
		{
			clover(&(GC[0]), &geo, &param, r, i, j, &M);
			trace_clover = retr(&M) + I*imtr(&M);
		
			trace_plaquettep = plaquettep_complex(&(GC[0]), &geo, &param, r, i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r, j), i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, r, i), i, j);
			trace_plaquettep += plaquettep_complex(&(GC[0]), &geo, &param, nnm(&geo, nnm(&geo, r, i), j), i, j);
		
			if(cabs(trace_clover-trace_plaquettep)<MIN_VALUE)
			{
				printf("Trace dir %d :	OK\n", i);
			}
			else
			{
				printf("	ERROR!!!!!!!!!!!\n\n	");
				printf("plane (%d, %d) clover = ", i, j);
				printfc(trace_clover);
				printf(", plaquettep = ");
				printfc(trace_plaquettep);
				printf("\n\n");
			}
		}
	}

	printf("\n\n");
	printf("VERIFY THAT TRANSLATIONS DON'T CHANGE THE MEAN PLAQUETTE ...\n\n");
	
	for(i=0; i<STDIM; i++)
	{
		plaquette(&(GC[0]), &geo, &param, &plaqs, &plaqt);			//mean plaquette
		conf_translation_dir(&(GC[0]), &geo, &param, i);			// translation
		plaquette(&(GC[0]), &geo, &param, &plaqs_new, &plaqt_new);	//mean plaquette new
		
		if(fabs(plaqs_new-plaqs)<MIN_VALUE)
		{
			printf("Space plaquette dir %d :	OK\n", i);
		}
		else
		{
			printf("	ERROR!!!!!!!!!!!\n\n	");
			printf("dir %d DeltaPlaqs = %g Plaqs = %g\n\n", i, plaqs_new-plaqs, plaqs);
		}
	
		if(fabs(plaqt_new-plaqt)<MIN_VALUE)
		{
			printf("Time plaquette dir %d :	OK\n", i);
		}
		else
		{
			printf("	ERROR!!!!!!!!!!!\n\n	");
			printf("dir %d DeltaPlaqt = %g Plaqt = %g\n\n", i, plaqt_new-plaqt, plaqt);
		}
	}
	
	printf("\n\n");
	printf("VERIFY THAT metropolis_single_swap SWAPS THE TWIST FACTORS ...\n\n");

	err=0;
	for(r=0; r<param.d_volume; r++)
	{
		for(j=0; j<param.d_n_planes; j++)
		{
			Z_a = GC[0].Z[r][j];
			Z_b = GC[1].Z[r][j];
			metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);
			
			if(cabs(GC[1].Z[r][j]-Z_a)>MIN_VALUE || cabs(GC[0].Z[r][j]-Z_b)> MIN_VALUE)
			{
				err = 1;
				printf("	ERROR!!!!!!!!!!!\n\n	");
				printf("Error swapping Z[%ld][%d]\n\n", r, j);
			}
		}
	}
	if(err==0) printf("Z[r][j] swaps :	OK\n");
	metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);
	
	printf("\n\n");
	printf("VERIFY THAT delta_action_swap IS CORRECT ...\n\n");
	
	err = 0;
	param.d_start = 1;
	init_gauge_conf_replica(&GC, &param);
	for(i=0; r<STDIM; i++)
		for(j=0; j<STDIM; j++)
			for(r=0; r<param.d_volume; r++)
			{
				K_a = (GC[0].C[r][i])*(GC[0].C[nnp(&geo, r, i)][j])*(GC[0].C[nnp(&geo, r, j)][i])*(GC[0].C[r][j]);
				K_b = (GC[1].C[r][i])*(GC[1].C[nnp(&geo, r, i)][j])*(GC[1].C[nnp(&geo, r, j)][i])*(GC[1].C[r][j]);
				re_tr_plaq_a = plaquettep(&(GC[0]),&geo,&param,r,i,j);
				re_tr_plaq_b = plaquettep(&(GC[1]),&geo,&param,r,i,j);
				action = -param.d_beta*(K_a*re_tr_plaq_a+K_b*re_tr_plaq_b);
				
				metropolis_single_swap(GC, 0, 1, 1.0, &acc_counters);
				
				K_a = (GC[0].C[r][i])*(GC[0].C[nnp(&geo, r, i)][j])*(GC[0].C[nnp(&geo, r, j)][i])*(GC[0].C[r][j]);
				K_b = (GC[1].C[r][i])*(GC[1].C[nnp(&geo, r, i)][j])*(GC[1].C[nnp(&geo, r, j)][i])*(GC[1].C[r][j]);
				re_tr_plaq_a = plaquettep(&(GC[0]),&geo,&param,r,i,j);
				re_tr_plaq_b = plaquettep(&(GC[1]),&geo,&param,r,i,j);
				swapped_action = -param.d_beta*(K_a*re_tr_plaq_a+K_b*re_tr_plaq_b);
				
				if(fabs(delta_action_swap(GC, &geo, &param, r, i, j, 0, 1)-(swapped_action-action))>MIN_VALUE )
				{
					err=1;
					printf("	ERROR!!!!!!!!!!!\n\n	");
					printf("Error site %ld, plane (%d,%d)\n\n", r, i, j);
				}
			}
	if(err==0) printf("delta_action_swap :	OK\n");

	free_replica(GC, &param);
	free_geometry(&geo, &param);
	free_rect_hierarc(most_update, clover_rectangle, &param);
	free_rect(&swap_rectangle);
	end_swap_acc_arrays(&acc_counters, &param);
	free_hierarc_params(&param);

	printf("\nTEST ENDED\n\n");

	return EXIT_SUCCESS;
	}

void conf_translation_dir(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param, int dir)
	{
	double aux;
	long s;
	Gauge_Conf aux_conf;

	// copy the conf in an auxiliary one (should be defined outside and passed to the function?), including the twist factors
	init_gauge_conf_from_gauge_conf(&aux_conf, GC, param); // now aux_conf=GC

	// translation in direction +dir, including the twist factors
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(s)
	#endif
	for(s=0;s<(param->d_n_planes)*(param->d_volume);s++)
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
	
void printfc(double complex z)
	{
	printf("%f%s%fi", creal(z), (cimag(z)>=0.0f)? "+":"", cimag(z));
	}
	
#endif
