/*
 * This is a sample-program for the vegas-algorithm. Compile the parallel
 * version using something similar to this (once you have the compiled
 * {p|n}vegas object-file):
 *   cc -o vegastest vegastest.c pvegas.o -lm -lpthread  (or equivalent) 
 * and the non-parallel (reference) version with:
 *   cc -o vegastest vegastest.c nvegas.o -lm
 * or the MPI-version equivalently (ask your administrator).  Please make 
 * sure to use the same compiler-switches as in production-runs or the 
 * test will be worthless!
 *
 * The test of the RNG is not really testing the RNG included in the
 * pvegas-sources because they all have different arguments for
 * technical reasons---instead it reemplements it.  It only checks
 * whether they are normal-distributed between 0 and 1.  (It is not the
 * algorithm which is being tested but what the compiler makes out of
 * it.)
 *
 * The test-function is a skewed gaussian that returns 1.0 when
 * integrated over the DIMENSION-dimensional hypercube even if WIDTH
 * gets close to 1.0 because the infinite real axis is mapped on the
 * finite interval (0.0, 1.0).  I do not claim that the choice of
 * vegas-calls is ideal, it is a good start---nothing more. Note
 * also, that the program is not at all suited for measuring timings
 * because testfun() is too lightweight. You need to supply your own
 * functions for that.  Good luck!
 */

#include <math.h>
#include <stdio.h>
#include "vegas.h"
#define DIMENSION 5
#define FUNCTIONS 4
#define RANDTEST 314159
#ifndef PI
#define PI     3.14159265358979323846
#endif

int t_gfsr_k;
unsigned int t_gfsr_m[SR_P];
double gfsr_norm;

#define WIDTH 0.2
void testfun(double x[DIMENSION], double f[FUNCTIONS])
{
  // double dummy, exponent, denominator;
  // int i;

  /* Calculate the principal integrand for which the grid will be optimized
     provided MAINFUNC has been defined as 0 in the Vegas modules: */
  // exponent = 0.0;
  // denominator = 1.0;
  // for (i=0; i<DIMENSION; ++i) {
  //   /* if you run into this one, you are really in trouble: */
  //   if (fabs(x[i])>=1) printf("Oh my god, x[%d]=%f\n",i,x[i]);
  //   /* this is a kluge for atanh() which doesn't exist on some systems: */
  //   dummy = 2*x[i]-1;
  //   dummy = 0.5*log((1+dummy)/(1-dummy));
  //   exponent -= dummy*dummy/2.0/WIDTH/WIDTH;
  //   dummy = 2*x[i]-1;
  //   denominator *= 1-dummy*dummy;
  // }
  // f[0] = exp(exponent)/denominator/pow(0.5*PI*WIDTH*WIDTH,(float)DIMENSION/2);
      f[0] = 1/(x[0]*x[0]+2*x[0]+2);
  /* It is slightly impractical that C does not have range-checking. Please
   * change FUNCTIONS to 2, 3 or 4 above and in the Vegas-module to integrate
   * additional functions along with the principal function f[0]. */
#if (FUNCTIONS>1)
  /* To first approximation this is a slight shift towards positive x[0]
   * where the peak is located at 0.5237 instead of at 0.5000 but leaving the
   * integral unchanged apart from overall constant exp(1): */
  f[1] = f[0] * 2.718281828 * (12.0*x[0]-5.0);
#endif
#if (FUNCTIONS>2)
  /* similar to acc1, but shifted towards positive x[1], normalized to pi: */
  f[2] = f[0] * 3.141592654 * (12.0*x[1]-5.0);
#endif
#if (FUNCTIONS>3)
  /* 0.0 in half the volume, f[0] in other half, integral normalized to 1.0: */
  f[3] = (x[2]<0.5) ? (f[0]*2.0) : (0.0);
#endif
  return;
}
double myrand(unsigned int w[], int *k)
{
  int j;

  (*k)++;
  if (*k >= SR_P) *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P) j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return((double)w[*k] * gfsr_norm);
}

void randtest()
{
  int i;
  unsigned int rdum = 42;
  double x, x0, sigma;

  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));
  for (i=0; i<SR_P; i++) {
    rdum = (1812433253*rdum + 314159265);
    t_gfsr_m[i] = rdum;
  }
  t_gfsr_k = -1;

  printf("Primitive test of the RNG:\n");
  printf(" expecting:  mean=0.500000   sigma=0.083333\n");
  sigma = x0 = 0;
  for (i=0; i<RANDTEST; i++) {
    x = myrand(t_gfsr_m,&t_gfsr_k);
    if (fabs(0.5-x)>0.5) printf("Oh my god, RN=%f\n",x);
    x0 += x;
  }
  printf(" calculated: mean=%f   ",x0/=RANDTEST);

  for (i=0; i<RANDTEST; i++) {
    x = myrand(t_gfsr_m,&t_gfsr_k);
    if (fabs(0.5-x)>0.5) printf("Oh my god, RN=%f\n",x);
    sigma += (x-x0)*(x-x0);
  }
  printf("sigma=%f   normalized with %g\n",sigma/=RANDTEST,gfsr_norm);
  return;
}

int main(int argc, char **argv)
{
  int i;
  double estim[FUNCTIONS];   /* estimators for integrals                     */
  double std_dev[FUNCTIONS]; /* standard deviations                          */
  double chi2a[FUNCTIONS];   /* chi^2/n                                      */
  double reg[2*DIMENSION];   /* integration domain                           */

  randtest();

#if defined(sun)    /* see README and set the argument below properly */
  if (thr_setconcurrency(8)) {
    perror("call of thr_setconcurrency() failed");
    exit(-1);
  };
#endif

  for (i=0; i<DIMENSION; i++) {
    reg[i] = 0.0;
    reg[i+DIMENSION] = 1.0;
  }

  /* set up the grid (init = 0) with 5 iterations of 1000 samples,
   * no need to compute additional accumulators (fcns = 1),
   * no parallelization yet (wrks = 1). */
  vegas(reg, DIMENSION, testfun,
        0, 1000, 5, NPRN_INPUT | NPRN_RESULT,
        1, 0, 1,
        estim, std_dev, chi2a);
  /* refine the grid (init = 1) with 5 iterations of 10000 samples,
   * collect in additional accumulators (fcns = FUNCTIONS),
   * two parallel workers (wrks = 2). */
  vegas(reg, DIMENSION, testfun,
        1, 10000, 5, NPRN_INPUT | NPRN_RESULT,
        FUNCTIONS, 0, 2,
        estim, std_dev, chi2a);
  /* final sample, inheriting previous results (init = 2) and
   * using 4 parallel workers (wrks = 4). */
  vegas(reg, DIMENSION, testfun,
        2, 100000, 2, NPRN_INPUT | NPRN_RESULT,
        FUNCTIONS, 0, 4,
        estim, std_dev, chi2a);

  printf ("Result: %g +/- %g\n", estim[0], std_dev[0]);
  // for (i=1; i<FUNCTIONS; ++i)
  //   printf("      ( %g +/- %g )\n", estim[i], std_dev[i]);
  /* Having set REPRO to 1 the numerical result should be (provided your
   * machine knows about IEEE)
   *      Principal integral  0.999142 +/- 0.000812909
   * 1st additional integral   2.71833 +/- 0.00406556 )
   * 2nd additional integral   3.13637 +/- 0.00468428 )
   * 3rd additional integral  0.998752 +/- 0.00115633 ) */
  return(0);
}
