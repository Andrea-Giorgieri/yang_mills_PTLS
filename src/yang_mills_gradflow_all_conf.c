#ifndef YM_GF_C
#define YM_GF_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    int stop;
	long step;
    double gftime;
	
    FILE *datafilep, *chiprimefilep, *topchar_tcorr_filep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file	
    readinput(in_file, &param);

    // this code has to start from saved conf.
    param.d_start=2;
	
	// not to overwrite files of runs with online gradient flow 
	strcpy(param.d_data_file, "_gradflow");
	strcpy(param.d_chiprime_file, "_gradflow");
	strcpy(param.d_topcharge_tcorr_file, "_gradflow");
	strcpy(param.d_log_file, "_gradflow");
	
    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &chiprimefilep, &topchar_tcorr_filep, &param);

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    time(&time1);
	if (param.d_saveconf_analysis_every == 0) stop = 1;
	else
	{
		stop = 0;
		step = ((int)(param.d_thermal/param.d_saveconf_analysis_every)+1)*param.d_saveconf_analysis_every;
	}
	
	while(stop == 0)
	{
		init_gauge_conf_step(&GC, &param, step, &stop);
		if (stop == 0)
		{
			perform_measures_localobs_with_gradflow(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tcorr_filep);
			step += param.d_saveconf_analysis_every;
		}
	}
    time(&time2);

    // close data file
    fclose(datafilep);
	if (param.d_chi_prime_meas==1) fclose(chiprimefilep);
	if (param.d_topcharge_tcorr_meas==1) fclose(topchar_tcorr_filep);

    // print simulation details
    print_parameters_gf(&param, time1, time2);

    // free gauge configurations
    free_gauge_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "#for gradient flow evolution\n");
    fprintf(fp, "gfstep      0.02    # integration step for gradient flow\n");
    fprintf(fp, "num_gfsteps 100     # number of integration steps for gradient flow\n");
	fprintf(fp, "gfstep_each 5       # compute observables every <gfstep_each> integration steps during the gradient flow\n");
	fprintf(fp, "\n");
	fprintf(fp, "#Simulation parameters\n");
	fprintf(fp, "thermal 0\n");
	fprintf(fp, "saveconf_analysis_every  100  # if 0 does not save, else save configurations for analysis every ... updates\n");
    fprintf(fp, "\n");
	fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
	fprintf(fp, "twist_file  twist.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
	fprintf(fp, "chiprime_data_file chi_prime_cool.dat\n");
	fprintf(fp, "topcharge_tcorr_file topo_tcorr_cool.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[500];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\tNum_levels (number of levels): %d\n", NLEVELS);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif

