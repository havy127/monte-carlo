/**
    Create by Ha Vy
    Date: 04/12/2019
    Topic: Monte Carlo 
 **/
#include <iostream>
#include <cmath>
#include "vegas.h"
using namespace std;

unsigned int gfsr_m[SR_P];           /* status-field for the GFSR-generator          */
int gfsr_k;                          /* pointer into status-field                    */
static int gfsr_not_initialized = 1; /* flag: need to initialize GFSR-field */
unsigned int rdum;                   /* linear congruential counter in kw_rand()     */
double gfsr_normn;                   /* will be set such that gfsr is normalized     */
int functions;   
/*
 * gfsr produces the random numbers once the starting values have been
 * initialized using gfsr_init. Origin from vegas
 */

/*
* Hàm tạo mãng ngẫu nhiên. 
*/
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
  return (rdum);
}

/**  Origin from vegas
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i, j;

  printf("Initializing SR-sequences with seed %ld\n", seed);
  gfsr_normn = (double)(1 / (pow(2.0, (long double)(8 * sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (j = 0; j < SR_P; j++)
    gfsr_m[j] = kw_rand(); /* initialize starting values */
  gfsr_k = -1;             /*         initialize pointer */
  gfsr_not_initialized = 0;
}