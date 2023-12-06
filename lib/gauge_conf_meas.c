#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/geometry.h"


#include<time.h> // DEBUG

// computation of the plaquette (1/NCOLOR the trace of, with twist factor) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j)
	{
	GAUGE_GROUP matrix;

	#ifdef DEBUG
	if(r >= param->d_volume)
	{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if(j >= STDIM || i >= STDIM)
	{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\nu\mu(x) = conj(Z_\mu\nu(x))
	
	return retr(&matrix);
	}


// computation of the plaquette (1/NCOLOR the trace of, with twist factor) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j)
	{
	GAUGE_GROUP matrix;
	
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if(j >= STDIM || i >= STDIM)
		{
		fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

	equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(&matrix, &(GC->lattice[r][i]));
	times_equal(&matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(&matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\mu\nu(x) = conj(Z_\nu\mu(x))
	
	return retr(&matrix)+I*imtr(&matrix);
	}


// computation of the plaquette with twist factor (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
	{
	#ifdef DEBUG
	if(r >= param->d_volume)
	{
		fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if(j >= STDIM || i >= STDIM)
	{
	fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}
	#else
	(void) param; // just to avoid warning at compile time
	#endif

//
//       ^ j
//       |   (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> i
//       r   (1)
//

	equal(matrix, &(GC->lattice[r][i]));
	times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
	times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
	times_equal_dag(matrix, &(GC->lattice[r][j]));
	
	//twist factor: Z_\mu\nu for clockwise plaquette with \mu < \nu, matrix is the anticlockwise plaquette
	times_equal_complex(matrix, GC->Z[r][dirs_to_si(j,i)]);	//Z_\mu\nu(x) = conj(Z_\nu\mu(x))
	}


// compute the four-leaf clover in position r, in the plane i,j and save it in M 
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux;
   long k, p;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= STDIM || j >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param;
   #endif

   zero(M);

//
//                   j ^
//                     |
//             (14)    |     (3)
//         +-----<-----++-----<-----+
//         |           ||           |
//         |           ||           |
//   (15)  V      (13) ^V (4)       ^ (2)
//         |           ||           |
//         |   (16)    || r   (1)   |
//    p    +----->-----++----->-----+------>   i
//         +-----<-----++-----<-----+
//         |    (9)    ||   (8)     |
//         |           ||           |
//    (10) V      (12) ^V (5)       ^ (7)
//         |           ||           |
//         |           ||           |
//         +------>----++----->-----+
//              (11)   k      (6)
//
   // avanti-avanti
   equal(&aux, &(GC->lattice[r][i]) );							// 1
   times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );		// 2
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );	// 3
   times_equal_dag(&aux, &(GC->lattice[r][j]) );				// 4
   times_equal_complex(&aux, GC->Z[r][dirs_to_si(i,j)]);		// twist anticlockwise
   plus_equal(M, &aux);

   k=nnm(geo, r, j);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][j]) );						// 5
   times_equal(&aux, &(GC->lattice[k][i]) );					// 6
   times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );		// 7
   times_equal_dag(&aux, &(GC->lattice[r][i]) );				// 8
   times_equal_complex(&aux, GC->Z[k][dirs_to_si(i,j)]);		// twist anticlockwise
   plus_equal(M, &aux);

   p=nnm(geo, r, i);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][i]) );								// 9
   times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );			// 10
   times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );				// 11
   times_equal(&aux, &(GC->lattice[k][j]) );							// 12
   times_equal_complex(&aux, GC->Z[nnm(geo, k, i)][dirs_to_si(i,j)]);	// twist anticlockwise
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][j]) );							// 13
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );	// 14
   times_equal_dag(&aux, &(GC->lattice[p][j]) );				// 15
   times_equal(&aux, &(GC->lattice[p][i]) );					// 16
   times_equal_complex(&aux, GC->Z[p][dirs_to_si(i,j)]);		// twist anticlockwise
   plus_equal(M, &aux);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, param, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, param, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy)
  {
  long r;
  double ris;

  ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<param->d_volume; r++)
     {
     int i, j;
     GAUGE_GROUP aux1, aux2;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           clover(GC, geo, param, r, i, j, &aux1);

           ta(&aux1);
           equal(&aux2, &aux1);
           times_equal(&aux1, &aux2);
           ris+=-NCOLOR*retr(&aux1)/16.0;
           }
        }
     }

  *energy=ris*param->d_inv_vol;
  }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly)
   {
   long rsp;
   double rep, imp;
   double complex tr;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }
      tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

      #if NCOLOR==1
        rep+=0.0;
      #else
        rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
      #endif

      imp+=0.0;
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }



// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_with_tracedef(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double *repoly,
                            double *impoly)
   {
   long rsp;
   double **rep, **imp;
   int j, err;
   long i;

   for(j=0;j<(int)floor(NCOLOR/2);j++)
      {
      repoly[j]=0.0;
      impoly[j]=0.0;
      }

   err=posix_memalign((void**)&rep, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&imp, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(i=0; i<param->d_space_vol; i++)
      {
      err=posix_memalign((void**)&(rep[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      err=posix_memalign((void**)&(imp[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      for(j=0; j<(int)floor(NCOLOR/2); j++)
         {
         rep[i][j] = 0.0;
         imp[i][j] = 0.0;
         }
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int k;
      GAUGE_GROUP matrix, matrix2;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(k=0; k<param->d_size[0]; k++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

       rep[rsp][0] = retr(&matrix);
       imp[rsp][0] = imtr(&matrix);

       equal(&matrix2, &matrix);

      for(k=1; k<(int)floor(NCOLOR/2.0); k++)
         {
         times_equal(&matrix2, &matrix);
         rep[rsp][k] = retr(&matrix2);
         imp[rsp][k] = imtr(&matrix2);
         }
      }

    for(j=0; j<(int)floor(NCOLOR/2); j++)
       {
       for(i=0; i<param->d_space_vol; i++)
          {
          repoly[j] += rep[i][j];
          impoly[j] += imp[i][j];
          }
       }

   for(j=0; j<(int)floor(NCOLOR/2.0); j++)
      {
      repoly[j] *= param->d_inv_space_vol;
      impoly[j] *= param->d_inv_space_vol;
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      free(rep[i]);
      free(imp[i]);
      }
   free(rep);
   free(imp);
   }


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					long r)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		(void) GC;
		(void) geo;
		(void) param;
		(void) r;
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	double ris;
	
	#if (STDIM==4 && NCOLOR>1)
	GAUGE_GROUP aux1, aux2, aux3;
	double real1, real2, loc_charge;
	const double chnorm=1.0/(128.0*PI*PI);
	int i, dir[4][3], sign;
	
	dir[0][0] = 0;
	dir[0][1] = 0;
	dir[0][2] = 0;
	
	dir[1][0] = 1;
	dir[1][1] = 2;
	dir[1][2] = 3;
	
	dir[2][0] = 2;
	dir[2][1] = 1;
	dir[2][2] = 1;
	
	dir[3][0] = 3;
	dir[3][1] = 3;
	dir[3][2] = 2;
	
	sign=-1;
	loc_charge=0.0;
	
	for(i=0; i<3; i++)
	{
		clover(GC, geo, param, r, dir[0][i], dir[1][i], &aux1);
		clover(GC, geo, param, r, dir[2][i], dir[3][i], &aux2);
		
		times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
		real1=retr(&aux3)*NCOLOR;
		
		times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
		real2=retr(&aux3)*NCOLOR;
		
		loc_charge+=((double) sign*(real1-real2));
		sign=-sign;
	}
	ris=(loc_charge*chnorm);
	#endif
	
	#if (STDIM==2 && NCOLOR==1)
	GAUGE_GROUP u1matrix;
	double angle;
	
	plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
	angle=atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;
	
	ris=angle;
	#endif
	
	return ris;
	}


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;
   long r;

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      ris+=loc_topcharge(GC, geo, param, r);
      }

   return ris;
   }

// sum loc_topcharge over STDIM-1 dirs and then the abs value of the result over the remaining dir
double topcharge_prime(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, int const dir)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	int i, err, cartcoord[STDIM];
	double ris, *tmp;
	long r;
	
	ris=0.0;
	err=posix_memalign((void**)&tmp, (size_t)DOUBLE_ALIGN, (size_t) param->d_size[dir] * sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for(i=0; i<param->d_size[dir]; i++) tmp[i] = 0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel num_threads(NTHREADS)
	#endif
	{
	int j;
	double *tmp_private;
	err=posix_memalign((void**)&tmp_private, (size_t)DOUBLE_ALIGN, (size_t) param->d_size[dir] * sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for(j=0; j<param->d_size[dir]; j++) tmp_private[j] = 0;
	
	#ifdef OPENMP_MODE
	#pragma omp for private(r, cartcoord) 
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		si_to_cart(cartcoord,r,param);
		tmp_private[cartcoord[dir]]+=loc_topcharge(GC, geo, param, r);
	}
	#ifdef OPENMP_MODE
	#pragma omp critical 
	#endif
	{
		for(j=0; j<param->d_size[dir]; j++) tmp[j] += tmp_private[j];
		free(tmp_private);
	}
	}

	for(i=0; i<param->d_size[dir]; i++) ris += fabs(tmp[i]);
	free(tmp);
	return ris;
	}

// chi^\prime = (1/8) int d^4x |x|^2 <q(x)q(0)> = < (1/8) int d^4x |x|^2 q(x) q(0) > = < G2 >
// This function computes the quantity (q(0)/8) sum_{x} d(x,0)^2 q(x) = a^2 G2, whose mean over the ensamble is a^2 chi^\prime
// d(x,y) = lattice distance between sites x and y keeping periodic boundary conditions into account (i.e., the shortest distance between x and y)
double topo_chi_prime(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	double ris=0.0, factor=0.125; // factor = 1/(2D) = 1/8
	long r;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris)
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		double d2 = square_distance(r, 0, param); // d(r,0)^2
		ris += d2 * loc_topcharge(GC, geo, param, r);
	}
	ris *=  loc_topcharge(GC, geo, param, 0) * factor; // ris *= q(0) / 8
	
	return ris;
}

void topcharge_timeslices(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, double *ris, int ncool, FILE *topchar_tcorr_filep)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	long r;
	int N_t = param->d_size[0];

	for (int i=0; i<N_t; i++) ris[i]=0.0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		int t = geo->d_timeslice[r];
		ris[t] += loc_topcharge(GC, geo, param, r);
	}

	fprintf(topchar_tcorr_filep, "%ld %d ", GC->update_index, ncool);
	for (int i=0; i<param->d_size[0]; i++) fprintf(topchar_tcorr_filep, " %.12g", ris[i]);
	fprintf(topchar_tcorr_filep, "\n");
}

void topcharge_timeslices_cooling(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, FILE *topchar_tcorr_filep)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	double *sum_q_timeslices;

	int err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}	

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf;
		int iter;

		// measure no cooling
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep); 
		// conf that will be cooled
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param); // helperconf is a copy of the configuration
		// measure with cooling
		for(iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, (iter+1)*param->d_coolsteps, topchar_tcorr_filep);
		}
		free_gauge_conf(&helperconf, param);
		fflush(topchar_tcorr_filep);
	}
	else // no cooling
	{
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
		fflush(topchar_tcorr_filep);
	}
	free(sum_q_timeslices);
	}

void topcharge_timeslices_gradflow(Gauge_Conf const * const GC,
									Geometry const * const geo,
									GParam const * const param,
									FILE *topchar_tcorr_filep)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	double *sum_q_timeslices;
	int err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	if(param->d_ngfsteps>0)	// if using gradient flow
	{	
		Gauge_Conf helperconf, help1, help2; 
		int count, meas_count;
		
		// measure no gradient flow
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, count, topchar_tcorr_filep);
			}
		}
		
		fflush(topchar_tcorr_filep);
		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
	}
	else	// no gradient flow
	{
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
		fflush(topchar_tcorr_filep);
	}
	free(sum_q_timeslices);
	}


// compute topological observables (Q, chi_prime) after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topo_obs_cooling(Gauge_Conf const * const GC,
					Geometry const * const geo,
					GParam const * const param,
					double *charge,
					double *chi_prime,
					double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

				ris=topo_chi_prime(&helperconf, geo, param);
				chi_prime[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, ris2, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
		 ris2=topo_chi_prime(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
				chi_prime[iter]=ris2;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }
	
void topo_obs_gradflow(Gauge_Conf const * const GC,
											Geometry const * const geo,
											GParam const * const param,
											double *charge,
											double *chi_prime,
											double *meanplaq)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_ngfsteps>0)	// if using gradient flow
	{	
		Gauge_Conf helperconf, help1, help2; 
		double plaqs, plaqt;
		int count, meas_count;
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				charge[meas_count]=topcharge(&helperconf, geo, param);
				chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
				plaquette(&helperconf, geo, param, &plaqs, &plaqt);
				#if(STDIM==4)
					meanplaq[meas_count]=0.5*(plaqs+plaqt);
				#else
					meanplaq[meas_count]=plaqt;
				#endif
			}
		}

		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
	}
	else	// no gradient flow
	{
		double plaqs, plaqt; 

		charge[0]=topcharge(GC, geo, param);
		chi_prime[0]=topo_chi_prime(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
		#if(STDIM==4)
			meanplaq[0]=0.5*(plaqs+plaqt);
		#else
			meanplaq[0]=plaqt;
		#endif
	}
	}
	
void topo_obs_clover_energy_gradflow(Gauge_Conf const * const GC,
											Geometry const * const geo,
											GParam const * const param,
											double *charge,
											double *chi_prime,
											double *clover_energy)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_ngfsteps>0)	// if using gradient flow
	{	
		Gauge_Conf helperconf, help1, help2;
		double tmp_energy;
		int count, meas_count;
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				charge[meas_count]=topcharge(&helperconf, geo, param);
				chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
				clover_disc_energy(&helperconf, geo, param, &tmp_energy);
				clover_energy[meas_count] = tmp_energy;
			}
		}

		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
	}
	else	// no gradient flow
	{
		double tmp_energy;
		
		charge[0]=topcharge(GC, geo, param);
		chi_prime[0]=topo_chi_prime(GC, geo, param);
		clover_disc_energy(GC, geo, param, &(clover_energy[0]));
	}
	}

/*---------------------------------------------*/
// OBSERVABLE NEEDED JUST TO CHECK HOW COOLING DESTROYS TOPOLOGICAL CORRELATIONS
void check_correlation_decay_cooling(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, double *ratio)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf; 
		double Q, satd;
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		for(int iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			Q = fabs(topcharge(&helperconf, geo, param));
			satd = sum_abs_topcharge_dens(&helperconf, geo, param);
			ratio[iter] = (satd-Q)/satd;
		}
		free_gauge_conf(&helperconf, param); 
	}
}

double sum_abs_topcharge_dens(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param)
{
	double sum=0.0;
	for (long r=0; r<(param->d_volume); r++)
	{
		sum += fabs(loc_topcharge(GC, geo, param, r));
	}
	return sum;
}

/*---------------------------------------------*/

// compute the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }

void topcharge_gradflow(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_ngfsteps>0)  // if using gradient flow
	{  
		Gauge_Conf helperconf, help1, help2; 
		double plaqs, plaqt;
		int count, meas_count;
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				charge[meas_count]=topcharge(&helperconf, geo, param);
				plaquette(&helperconf, geo, param, &plaqs, &plaqt);
				#if(STDIM==4)
					meanplaq[meas_count]=0.5*(plaqs+plaqt);
				#else
					meanplaq[meas_count]=plaqt;
				#endif

			}
		}

		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
	}
	else   // no gradient flow
	{
		double plaqs, plaqt; 
		int iter;

		charge[0]=topcharge(GC, geo, param);
		plaquette(GC, geo, param, &plaqs, &plaqt);
		#if(STDIM==4)
			meanplaq[0]=0.5*(plaqs+plaqt);
		#else
			meanplaq[0]=plaqt;
		#endif
	}
	}
	
void topcharge_clover_energy_gradflow(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *clover_energy)
	{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_ngfsteps>0)  // if using gradient flow
	{  
		Gauge_Conf helperconf, help1, help2;
		double tmp_energy;
		int count, meas_count;
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				charge[meas_count]=topcharge(&helperconf, geo, param);
				clover_disc_energy(&helperconf, geo, param, &tmp_energy);
				clover_energy[meas_count] = tmp_energy;
			}
		}

		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
	}
	else   // no gradient flow
	{
		double tmp_energy;
		
		charge[0]=topcharge(GC, geo, param);
		clover_disc_energy(GC, geo, param, &tmp_energy);
		clover_energy[0] = tmp_energy;
	}
	}


// compute the correlator of the local topological charge
// after "ncool" cooling steps up to spatial distance "dist"
void loc_topcharge_corr(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    int ncool,
                    int dist,
                    double *ris)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double *topch;
   long r;
   int err, i;

   err=posix_memalign((void**) &(topch), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating memory! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   // compute the local topological charge
   if(ncool>0)
     {
     Gauge_Conf helperconf;

     // helperconf is a copy of GC
     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);

     // cool helperconf
     cooling(&helperconf, geo, param, ncool);

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(&helperconf, geo, param, r);
        }

     // free helperconf
     free_gauge_conf(&helperconf, param);
     }
   else
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(GC, geo, param, r);
        }
     }

   // compute correlators
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(i)
   #endif
   for(i=0; i<dist; i++)
      {
      double aux;
      long r1, r2;
      int j, dir;

      ris[i]=0.0;

      for(r1=0; r1<param->d_volume; r1++)
         {
         aux=0.0;

         for(dir=1; dir<STDIM; dir++)
            {
            r2=r1;
            for(j=0; j<i; j++)
               {
               r2=nnp(geo, r2, dir);
               }
            aux+=topch[r2];
            }
         aux/=(double)(STDIM-1);

         ris[i]+=aux*topch[r1];
         }
      ris[i]*=param->d_inv_vol;
      }

   // free memory
   free(topch);
   }


void perform_measures_localobs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double plaqs, plaqt, polyre, polyim, *charge, *chi_prime, *meanplaq, charge_nocooling, chi_prime_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);
		 
     charge_nocooling=topcharge(GC, geo, param);
		 if (param->d_chi_prime_meas == 1 ) chi_prime_nocooling=topo_chi_prime(GC, geo, param);

		 // refresh topological charge of periodic replica (only for multicanonic)
		 GC->stored_topo_charge = charge_nocooling;

     fprintf(datafilep, "%ld %.12g %.12g %.12g %.12g %.12g ", GC->update_index, plaqs, plaqt, polyre, polyim, charge_nocooling);
		 if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_chi_prime_meas == 1)
		 {
     	err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     	if(err!=0)
       	{
       	fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       	exit(EXIT_FAILURE);
       	}
		 }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_topcharge_tcorr_meas == 1 ) topcharge_timeslices_cooling(GC, geo, param, topchar_tcorr_filep);
		 else {(void) topchar_tcorr_filep;}
     if (param->d_chi_prime_meas == 1 ) topo_obs_cooling(GC, geo, param, charge, chi_prime, meanplaq);
		 else topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
				if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_coolsteps, chi_prime[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);
		 if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);

     free(charge);
		 if (param->d_chi_prime_meas == 1 ) free(chi_prime);
     else 
			{
				(void) chiprimefilep;
			}
		 free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);
   #endif
   }
	
void perform_measures_localobs_with_gradflow(Gauge_Conf *GC,
											Geometry const * const geo,
											GParam const * const param,
											FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int i, j, err, gradflowrepeat;
	double plaqs, plaqt, polyre, polyim, clover_energy_nogradflow, charge_nogradflow, chi_prime_nogradflow, charge_prime_nogradflow[STDIM];
	double *meanplaq, *clover_energy, *charge, *sum_q_timeslices, *chi_prime, **charge_prime;
	
	// meas no gradflow
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) charge_nogradflow=topcharge(GC, geo, param);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	if (param->d_chi_prime_meas == 1 ) chi_prime_nogradflow=topo_chi_prime(GC, geo, param);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) charge_prime_nogradflow[i]=topcharge_prime(GC, geo, param, i);
	if (param->d_topcharge_tcorr_meas == 1 )
	{
		err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
	}
	else {(void) topchar_tcorr_filep;}
	
	// refresh topological charge of periodic replica (only for multicanonic)
	GC->stored_topo_charge = charge_nogradflow;

	// print meas no gradflow (topcharge_tcorr_timeslices already printed by topcharge_timeslices())
	fprintf(datafilep, "%ld ", GC->update_index);
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge_nogradflow);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nogradflow);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) fprintf(datafilep, "%.12g ", charge_prime_nogradflow[i]);
	
	// meas gradflow
	gradflowrepeat = (int)(param->d_ngfsteps/param->d_gf_meas_each);
	if (gradflowrepeat > 0)
	{
		Gauge_Conf helperconf, help1, help2;
		int count, meas_count;
		
		if (param->d_plaquette_meas == 1)
		{
			err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_clover_energy_meas == 1)
		{
			err=posix_memalign((void**)&clover_energy, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_charge_meas == 1)
		{
			err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
				
		if (param->d_chi_prime_meas == 1)
		{
			err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		
		if (param->d_charge_prime_meas == 1)
		{
			err=posix_memalign((void**)&charge_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double*));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			for(i=0; i<gradflowrepeat; i++)
			{
				err=posix_memalign((void**)&(charge_prime[i]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
				if(err!=0)
				{
					fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				if (param->d_plaquette_meas == 1 ) 
				{
					plaquette(&helperconf, geo, param, &plaqs, &plaqt);
					#if(STDIM==4)
						meanplaq[meas_count]=0.5*(plaqs+plaqt);
					#else
						meanplaq[meas_count]=plaqt;
					#endif
				}
				if (param->d_clover_energy_meas == 1 ) clover_disc_energy(&helperconf, geo, param, &clover_energy[meas_count]);
				if (param->d_charge_meas == 1 ) charge[meas_count]=topcharge(&helperconf, geo, param);
				if (param->d_topcharge_tcorr_meas == 1 )
				{
					topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, count, topchar_tcorr_filep);
				}				
				if (param->d_chi_prime_meas == 1) chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
				if (param->d_charge_prime_meas == 1) for (i=0; i<STDIM; i++) charge_prime[meas_count][i]=topcharge_prime(&helperconf, geo, param, i);
			}
		}
		
		// print meas gradflow
		for(i=0; i<gradflowrepeat; i++)
		{
			if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g ", meanplaq[i]);
			if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy[i]);
			if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge[i]);
			if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_ngfsteps, chi_prime[i]);
			if (param->d_charge_prime_meas == 1 ) for (j=0; j<STDIM; j++) fprintf(datafilep, "%.12g ", charge_prime[i][j]);
		}
		
		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
		if (param->d_plaquette_meas == 1 ) free(meanplaq);
		if (param->d_clover_energy_meas == 1 ) free(clover_energy);
		if (param->d_charge_meas == 1 ) free(charge);
		if (param->d_topcharge_tcorr_meas == 1 ) free(sum_q_timeslices);
		if (param->d_chi_prime_meas == 1 ) free(chi_prime);
		if (param->d_charge_prime_meas == 1 )
		{
			for(i=0; i<gradflowrepeat; i++) free(charge_prime[i]);
			free(charge_prime);
		}
	}
	
	fprintf(datafilep, "\n");
	fflush(datafilep);
	if (param->d_topcharge_tcorr_meas == 1 ) fflush(topchar_tcorr_filep);
	if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);
	
	#else

	double plaqs, plaqt, clover_energy, polyre, polyim;
	
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	fprintf(datafilep, "\n");
	fflush(datafilep);
	#endif
	}

void perform_measures_localobs_with_adaptive_gradflow(Gauge_Conf *GC,
											Geometry const * const geo,
											GParam const * const param,
											FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int i, j, err, gradflowrepeat;
	double plaqs, plaqt, polyre, polyim, clover_energy_nogradflow, charge_nogradflow, chi_prime_nogradflow, charge_prime_nogradflow[STDIM];
	double *meanplaq, *clover_energy, *charge, *sum_q_timeslices, *chi_prime, **charge_prime;
	
	// meas no gradflow
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) charge_nogradflow=topcharge(GC, geo, param);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	if (param->d_chi_prime_meas == 1 ) chi_prime_nogradflow=topo_chi_prime(GC, geo, param);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) charge_prime_nogradflow[i]=topcharge_prime(GC, geo, param, i);
	if (param->d_topcharge_tcorr_meas == 1 )
	{
		err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
	}
	else {(void) topchar_tcorr_filep;}
	
	// refresh topological charge of periodic replica (only for multicanonic)
	GC->stored_topo_charge = charge_nogradflow;

	// print meas no gradflow (topcharge_tcorr_timeslices already printed by topcharge_timeslices())
	fprintf(datafilep, "%ld ", GC->update_index);
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge_nogradflow);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nogradflow);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) fprintf(datafilep, "%.12g ", charge_prime_nogradflow[i]);
	
	// meas gradflow
	gradflowrepeat = (int)(param->d_agf_length/param->d_agf_meas_each);
	if (gradflowrepeat > 0)
		{
		Gauge_Conf helperconf_old, helperconf, help1, help2, help3;
		int meas_count, accepted;
		double gftime_old, gftime, gftime_step;
		
		if (param->d_plaquette_meas == 1)
		{
			err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_clover_energy_meas == 1)
		{
			err=posix_memalign((void**)&clover_energy, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_charge_meas == 1)
		{
			err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
				
		if (param->d_chi_prime_meas == 1)
		{
			err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		
		if (param->d_charge_prime_meas == 1)
		{
			err=posix_memalign((void**)&charge_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double*));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			for(i=0; i<gradflowrepeat; i++)
			{
				err=posix_memalign((void**)&(charge_prime[i]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
				if(err!=0)
				{
					fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		
		init_gauge_conf_from_gauge_conf(&helperconf_old, GC, param);
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);
		init_gauge_conf_from_gauge_conf(&help3, GC, param);

		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		meas_count = 0;
		while(meas_count < gradflowrepeat)
			{
			gradflow_RKstep_adaptive(&helperconf, &helperconf_old, &help1, &help2, &help3, geo, param, &gftime, &gftime_step, &accepted);
			if (accepted == 1 && fabs(gftime - param->d_agf_meas_each*(meas_count+1)) < param->d_agf_time_bin + MIN_VALUE) 	//step accepted, perform measures
				{
				if (param->d_plaquette_meas == 1 ) 
					{
					plaquette(&helperconf, geo, param, &plaqs, &plaqt);
					#if(STDIM==4)
					meanplaq[meas_count]=0.5*(plaqs+plaqt);
					#else
					meanplaq[meas_count]=plaqt;
					#endif
					}
				if (param->d_clover_energy_meas == 1 ) clover_disc_energy(&helperconf, geo, param, &clover_energy[meas_count]);
				if (param->d_charge_meas == 1 ) charge[meas_count]=topcharge(&helperconf, geo, param);
				if (param->d_topcharge_tcorr_meas == 1 )
					{
					topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, meas_count, topchar_tcorr_filep);
					}				
				if (param->d_chi_prime_meas == 1) chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
				if (param->d_charge_prime_meas == 1) for (i=0; i<STDIM; i++) charge_prime[meas_count][i]=topcharge_prime(&helperconf, geo, param, i);
				meas_count = meas_count + 1;
				}
			if ((gftime + gftime_step - param->d_agf_meas_each*(meas_count+1)) > param->d_agf_time_bin ) //adapt step to the time of next measure
				{
				gftime_step = param->d_agf_meas_each*(meas_count+1) - gftime;
				}
			}
		
		// print meas gradflow
		for(i=0; i<gradflowrepeat; i++)
			{
			if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g ", meanplaq[i]);
			if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy[i]);
			if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge[i]);
			if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_ngfsteps, chi_prime[i]);
			if (param->d_charge_prime_meas == 1 ) for (j=0; j<STDIM; j++) fprintf(datafilep, "%.12g ", charge_prime[i][j]);
			}
		
		free_gauge_conf(&helperconf_old, param);
		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
		free_gauge_conf(&help3, param);
		if (param->d_plaquette_meas == 1 ) free(meanplaq);
		if (param->d_clover_energy_meas == 1 ) free(clover_energy);
		if (param->d_charge_meas == 1 ) free(charge);
		if (param->d_topcharge_tcorr_meas == 1 ) free(sum_q_timeslices);
		if (param->d_chi_prime_meas == 1 ) free(chi_prime);
		if (param->d_charge_prime_meas == 1 )
			{
			for(i=0; i<gradflowrepeat; i++) free(charge_prime[i]);
			free(charge_prime);
			}
		}
	
	fprintf(datafilep, "\n");
	fflush(datafilep);
	if (param->d_topcharge_tcorr_meas == 1 ) fflush(topchar_tcorr_filep);
	if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);
	
	#else

	double plaqs, plaqt, clover_energy, polyre, polyim;
	
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	fprintf(datafilep, "\n");
	fflush(datafilep);
	#endif
	}

void perform_measures_localobs_with_adaptive_gradflow_debug(Gauge_Conf *GC,
											Geometry const * const geo,
											GParam const * const param,
											FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep, FILE *step_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int i, j, err, gradflowrepeat;
	double plaqs, plaqt, polyre, polyim, clover_energy_nogradflow, charge_nogradflow, chi_prime_nogradflow, charge_prime_nogradflow[STDIM];
	double *meanplaq, *clover_energy, *charge, *sum_q_timeslices, *chi_prime, **charge_prime;
	
	// meas no gradflow
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) charge_nogradflow=topcharge(GC, geo, param);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	if (param->d_chi_prime_meas == 1 ) chi_prime_nogradflow=topo_chi_prime(GC, geo, param);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) charge_prime_nogradflow[i]=topcharge_prime(GC, geo, param, i);
	if (param->d_topcharge_tcorr_meas == 1 )
	{
		err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
	}
	else {(void) topchar_tcorr_filep;}
	
	// refresh topological charge of periodic replica (only for multicanonic)
	GC->stored_topo_charge = charge_nogradflow;

	// print meas no gradflow (topcharge_tcorr_timeslices already printed by topcharge_timeslices())
	fprintf(datafilep, "%ld ", GC->update_index);
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy_nogradflow);
	if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge_nogradflow);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nogradflow);
	if (param->d_charge_prime_meas == 1 ) for (i=0; i<STDIM; i++) fprintf(datafilep, "%.12g ", charge_prime_nogradflow[i]);
	
	// meas gradflow
	gradflowrepeat = (int)(param->d_agf_length/param->d_agf_meas_each);
	if (gradflowrepeat > 0)
		{
		Gauge_Conf helperconf_old, helperconf, help1, help2, help3;
		int meas_count, accepted;
		double gftime_old, gftime, gftime_step;
		
		if (param->d_plaquette_meas == 1)
		{
			err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_clover_energy_meas == 1)
		{
			err=posix_memalign((void**)&clover_energy, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
			
		if (param->d_charge_meas == 1)
		{
			err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
				
		if (param->d_chi_prime_meas == 1)
		{
			err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		
		if (param->d_charge_prime_meas == 1)
		{
			err=posix_memalign((void**)&charge_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double*));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			for(i=0; i<gradflowrepeat; i++)
			{
				err=posix_memalign((void**)&(charge_prime[i]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
				if(err!=0)
				{
					fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		
		init_gauge_conf_from_gauge_conf(&helperconf_old, GC, param);
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);
		init_gauge_conf_from_gauge_conf(&help3, GC, param);

		// gradflow starts
		gftime = 0.0;
		gftime_step = param->d_agf_step;
		meas_count = 0;
		while(meas_count < gradflowrepeat)
			{
			fprintf(step_filep, "%ld %.12g %.12g ", GC->update_index, gftime, gftime_step);
			gradflow_RKstep_adaptive_debug(&helperconf, &helperconf_old, &help1, &help2, &help3, geo, param, &gftime, &gftime_step, &accepted, step_filep);
			if (accepted == 1 && fabs(gftime - param->d_agf_meas_each*(meas_count+1)) <= param->d_agf_time_bin ) 	//step accepted, perform measures
				{
				if (param->d_plaquette_meas == 1 ) 
					{
					plaquette(&helperconf, geo, param, &plaqs, &plaqt);
					#if(STDIM==4)
					meanplaq[meas_count]=0.5*(plaqs+plaqt);
					#else
					meanplaq[meas_count]=plaqt;
					#endif
					}
				if (param->d_clover_energy_meas == 1 ) clover_disc_energy(&helperconf, geo, param, &clover_energy[meas_count]);
				if (param->d_charge_meas == 1 ) charge[meas_count]=topcharge(&helperconf, geo, param);
				if (param->d_topcharge_tcorr_meas == 1 )
					{
					topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, meas_count, topchar_tcorr_filep);
					}				
				if (param->d_chi_prime_meas == 1) chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
				if (param->d_charge_prime_meas == 1) for (i=0; i<STDIM; i++) charge_prime[meas_count][i]=topcharge_prime(&helperconf, geo, param, i);
				meas_count = meas_count + 1;
				}
			if ((gftime + gftime_step - param->d_agf_meas_each*(meas_count+1)) > param->d_agf_time_bin ) //adapt step to the time of next measure
				{
				gftime_step = param->d_agf_meas_each*(meas_count+1) - gftime;
				}
			fprintf(step_filep, "\n");
			}
		
		// print meas gradflow
		for(i=0; i<gradflowrepeat; i++)
			{
			if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g ", meanplaq[i]);
			if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy[i]);
			if (param->d_charge_meas == 1 ) fprintf(datafilep, "%.12g ", charge[i]);
			if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_ngfsteps, chi_prime[i]);
			if (param->d_charge_prime_meas == 1 ) for (j=0; j<STDIM; j++) fprintf(datafilep, "%.12g ", charge_prime[i][j]);
			}
		
		free_gauge_conf(&helperconf_old, param);
		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
		free_gauge_conf(&help3, param);
		if (param->d_plaquette_meas == 1 ) free(meanplaq);
		if (param->d_clover_energy_meas == 1 ) free(clover_energy);
		if (param->d_charge_meas == 1 ) free(charge);
		if (param->d_topcharge_tcorr_meas == 1 ) free(sum_q_timeslices);
		if (param->d_chi_prime_meas == 1 ) free(chi_prime);
		if (param->d_charge_prime_meas == 1 )
			{
			for(i=0; i<gradflowrepeat; i++) free(charge_prime[i]);
			free(charge_prime);
			}
		}
	
	fprintf(datafilep, "\n");
	fflush(datafilep);
	fflush(step_filep);
	if (param->d_topcharge_tcorr_meas == 1 ) fflush(topchar_tcorr_filep);
	if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);
	
	#else

	double plaqs, plaqt, clover_energy, polyre, polyim;
	
	if (param->d_plaquette_meas == 1 ) plaquette(GC, geo, param, &plaqs, &plaqt);
	if (param->d_clover_energy_meas == 1 ) clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	
	if (param->d_plaquette_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
	if (param->d_clover_energy_meas == 1 ) fprintf(datafilep, "%.12g ", clover_energy);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	fprintf(datafilep, "\n");
	fflush(datafilep);
	#endif
	}
	
void perform_measures_localobs_clover_energy_with_gradflow(Gauge_Conf *GC,
											Geometry const * const geo,
											GParam const * const param,
											FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep)
	{
	#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
	int i, err, gradflowrepeat;
	double clover_energy_nogradflow, polyre, polyim, *charge, *sum_q_timeslices, *chi_prime, *clover_energy, charge_nogradflow, chi_prime_nogradflow;
	
	// meas no gradflow
	clover_disc_energy(GC, geo, param, &clover_energy_nogradflow);
	charge_nogradflow=topcharge(GC, geo, param);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	if (param->d_chi_prime_meas == 1 ) chi_prime_nogradflow=topo_chi_prime(GC, geo, param);
	if (param->d_topcharge_tcorr_meas == 1 )
	{
		err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
	}
	else {(void) topchar_tcorr_filep;}
	
	// refresh topological charge of periodic replica (only for multicanonic)
	GC->stored_topo_charge = charge_nogradflow;

	// print meas no gradflow
	fprintf(datafilep, "%ld %.12g %.12g ", GC->update_index, clover_energy_nogradflow, charge_nogradflow);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nogradflow);
	
	// meas gradflow
	gradflowrepeat = (int)(param->d_ngfsteps/param->d_gf_meas_each);
	if (gradflowrepeat > 0)
	{
		Gauge_Conf helperconf, help1, help2;
		int count, meas_count;
		
		err=posix_memalign((void**)&clover_energy, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		
		err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		
		if (param->d_chi_prime_meas == 1)
		{
			err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) gradflowrepeat * sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		init_gauge_conf_from_gauge_conf(&help1, GC, param);
		init_gauge_conf_from_gauge_conf(&help2, GC, param);

		// count starts from 1 to avoid problems with %
		for(count=1; count < (param->d_ngfsteps+1); count++)
		{
			gradflow_RKstep(&helperconf, &help1, &help2, geo, param, param->d_gfstep);
			
			if ( (count % param->d_gf_meas_each) == 0)
			{
				meas_count = count/param->d_gf_meas_each-1;
				clover_disc_energy(&helperconf, geo, param, &clover_energy[meas_count]);
				charge[meas_count]=topcharge(&helperconf, geo, param);
				if (param->d_topcharge_tcorr_meas == 1 )
				{
					topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, count, topchar_tcorr_filep);
				}				
				if (param->d_chi_prime_meas == 1) chi_prime[meas_count]=topo_chi_prime(&helperconf, geo, param);
			}
		}
		
		// print meas gradflow
		for(i=0; i<gradflowrepeat; i++)
		{
			fprintf(datafilep, "%.12g %.12g ", clover_energy[i], charge[i]);
			if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_ngfsteps, chi_prime[i]);
		}
		
		free_gauge_conf(&helperconf, param);
		free_gauge_conf(&help1, param);
		free_gauge_conf(&help2, param);
		free(clover_energy);
		free(charge);
		if (param->d_topcharge_tcorr_meas == 1 ) free(sum_q_timeslices);
		if (param->d_chi_prime_meas == 1 ) free(chi_prime);
	}
	
	fprintf(datafilep, "\n");
	fflush(datafilep);
	if (param->d_topcharge_tcorr_meas == 1 ) fflush(topchar_tcorr_filep);
	if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);
	
	#else

	double clover_energy, polyre, polyim;
	
	clover_disc_energy(GC, geo, param, &clover_energy);
	if (param->d_polyakov_meas == 1 ) polyakov(GC, geo, param, &polyre, &polyim);
	
	fprintf(datafilep, "%.12g ", clover_energy);
	if (param->d_polyakov_meas == 1 ) fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
	fprintf(datafilep, "\n");
	fflush(datafilep);
	#endif
	}

// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )

     int i, err;
     double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1
     double *charge, *meanplaq, charge_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov_with_tracedef(GC, geo, param, polyre, polyim);
     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

     for(i=0; i<(int)floor(NCOLOR/2); i++)
        {
        fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
        }
     fprintf(datafilep, "%.12g ", charge_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);

     free(charge);
     free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);

   #endif
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double complex poly_corr;
  double poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }
        poly_array[r]=retr(&matrix)+I*imtr(&matrix);
        }

     // average correlator computation
     poly_corr=0.0+I*0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr += poly_array[r]*conj(poly_array[r2]);
        poly_corr_abs += cabs(poly_array[r]*conj(poly_array[r2]));
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_corr_fluct+=cabs( poly_array[r]*conj(poly_array[r2]) - poly_corr );
        }
     poly_corr_fluct*=param->d_inv_space_vol;


     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   int i, err;
   long r;
   double complex poly_corr;
   double poly_corr_abs, poly_corr_fluct;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorr(GC,
                       geo,
                       param,
                       param->d_size[0]);
   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct += cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*= sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g ", poly_corr_abs);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorr(GC,
                geo,
                param,
                param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorr(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorr(GC, geo, param, datafilep);
   #endif
   }


// to optimize the number of hits to be used in multilevel for the adjoint representation
void optimize_multihit_polycorradj(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double poly_corr, poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex tr;
  double *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop in the adjoint representation computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }

        //trace of the matix in the fundamental representation
        tr=NCOLOR*(retr(&matrix)+I*imtr(&matrix));

        //trace of the matrix in adjoint representation
        poly_array[r]=cabs(tr*conj(tr))-1.0;
        #if NCOLOR != 1
          poly_array[r]/=(NCOLOR*NCOLOR-1);
        #endif
        }

     // average correlator computation
     poly_corr=0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr+= poly_array[r]*poly_array[r2];
        poly_corr_abs+=fabs(poly_array[r]*poly_array[r2]);
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr_fluct+=fabs(poly_array[r]*poly_array[r2]-poly_corr);
        }
     poly_corr_fluct*=param->d_inv_space_vol;

     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel (adjoint representation)
void optimize_multilevel_polycorradj(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr, poly_corr_abs, poly_corr_fluct;
   double *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorradj(GC,
                          geo,
                          param,
                          param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=fabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=fabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator in the adjoint representation with the multilevel algorithm
void perform_measures_polycorradj(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorradj(GC,
                            geo,
                            param,
                            param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorradj(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorradj(GC, geo, param, datafilep);
   #endif
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr_abs, poly_corr_fluct;
   double complex poly_corr;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // average
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuation
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalization
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_ml_level0_repeat);
   poly_corr_fluct*=sqrt(param->d_ml_level0_repeat);

   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *GC,
                                    GParam const * const param,
                                    FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
      optimize_multilevel_polycorr_long(GC, param, datafilep);
   #else
     double ris;
     long r;
     int i;

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   }


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_disc(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_conn(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif
