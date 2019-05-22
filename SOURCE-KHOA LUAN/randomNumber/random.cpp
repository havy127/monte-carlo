#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "radmom.h"
using namespace std;



static int gfsr_not_initialized = 1; /* flag: need to initialize GFSR-field */
int gfsr_k;                          /* pointer into status-field                    */
unsigned int gfsr_m[SR_P];           /* status-field for the GFSR-generator          */
unsigned int rdum;                   /* linear congruential counter in kw_rand()     */
double gfsr_normn;                   /* will be set such that gfsr is normalized     */

unsigned int kw_rand(void)
{
  rdum = (1812433253 * rdum + 314159265);
  cout << "\n rdum:  " << rdum;
  return (rdum);
}

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
  return ((double)w[*k] * gfsr_normn );
}
void randomNumberTest()
{
  
   double num;
    gfsr_init((long)time(NULL));
    num = gfsr_rand(gfsr_m, &gfsr_k);
   cout << "\n" <<rdum << endl; 

}
