/*
 * NAME
 *   eduvegas.c
 *   Implementation of G.P.Lepage's VEGAS-algorithm and
 *      Richard Kreckel, ThEP, Univ. Mainz, October 1996 - October 2000 
 *
 *
 * AIM eduvegas.c
 *   modified for education purposes: Explain statistics knowlegde in vegas
 *
 * STATUS
 *   Under modifing .... current time: November 2017,
 *
 * SYNOPSIS
 *   void vegas(double regn[], int ndim, void (*fxn)(double x[], double f),
 *              int init, unsigned long ncall, int itmx, int nprn,
 *              double tgral[], double sd[], double chi2a[]);
 *
 *     regn[]: array specifying the region to be integrated, 2*ndim entries
 *     ndim: dimensionality of space
 *     (*fxn)(x[],f[]): pointer to function to be evaluated (must be MT-safe!)
 *     init: initialization level (start with 0, then 1, later 2)
 *     ncall: number of samples points per iteration
 *     itmx: number of iterations in one call
 *     nprn: bit field, see constants NPRN_* below
 *     tgral[]: pointer to estimate of result (maybe array)
 *     sd[]: pointer to estimate of standard deviation (maybe array)
 *     chi2a[]: pointer to chi-squared over ndim (maybe array)*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vegas.h"

#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define REPRO 1      /* 0 = default, others used for comparison      */

unsigned int gfsr_m[SR_P];           /* status-field for the GFSR-generator          */
int gfsr_k;                          /* pointer into status-field                    */
static int gfsr_not_initialized = 1; /* flag: need to initialize GFSR-field */
unsigned int rdum;                   /* linear congruential counter in kw_rand()     */
double gfsr_normn;                   /* will be set such that gfsr is normalized     */
int functions;                       /* copy of (*ctl).fcns                          */



double gfsr_rand(unsigned int w[], int *k)
{
  int j;

  (*k)++;
  if (*k >= SR_P)
    *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P)
    j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return ((double)w[*k] * gfsr_normn);
}

/**  Origin from vegas
 * Simple linear congruential generator used to initialize SR.
 * The multiplier and increment were provided by KOMA, Univ.Mainz.
 * (kw stands for Kalos & Whitlock)
 */
/*
* Giã số ngẫu nhiên
*/
unsigned int kw_rand(void)
{
  rdum = (1812433253 * rdum + 314159265);
  printf ("\n rdum: %d",rdum);
  return (rdum);
}

/**  Origin from vegas
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i, j,k;

  printf("Initializing SR-sequences with seed %ld\n", seed);
  gfsr_normn = (double)(1 / (pow(2.0, (long double)(8 * sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (j = 0; j < SR_P; j++)
  gfsr_m[j] = kw_rand(); /* initialize starting values */
  gfsr_k = -1;             /*         initialize pointer */
  gfsr_not_initialized = 0;
  // for (k = 0; k < SR_P; k++)
  // printf ("\n gfsr_m[%d]: %d",k,gfsr_m[k]);
}
//=======================================================


//=======================================================

void randomNumberTest()
{
  
   double num;
    gfsr_init((long)time(NULL));
    num = gfsr_rand(gfsr_m, &gfsr_k);

   printf ("\nRandom number: %lf",num); 

}
//=======================================================


#undef SR_P
#undef SR_Q
#undef REPRO
#undef TINY
