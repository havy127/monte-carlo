/**
 * Created by PhpStorm.
 * User: Huynh Thi Ha Vy
 * Date: 03/04/2019
 * Time: 4:24 PM
 */

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
using namespace std;
#include "vegas.h"

#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define REPRO 1      /* 0 = default, others used for comparison      */

unsigned int gfsr_m[SR_P];           /* status-field for the GFSR-generator          */
int gfsr_k;                          /* pointer into status-field                    */
static int gfsr_not_initialized = 1; /* flag: need to initialize GFSR-field */
unsigned int rdum;                   /* linear congruential counter in kw_rand()     */
double gfsr_normn;                   /* will be set such that gfsr is normalized     */
int functions;  


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
  
  cout << "Initializing SR-sequences with seed" << endl << seed <<endl;
  gfsr_normn = (double)(1 / (pow(2.0, (long double)(8 * sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (j = 0; j < SR_P; j++)
    gfsr_m[j] = kw_rand(); /* initialize starting values */
  gfsr_k = -1;             /*         initialize pointer */
  gfsr_not_initialized = 0;
}