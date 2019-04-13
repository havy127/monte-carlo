#ifndef _VEGASCLASS_
#define _VEGASCLASS_

#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define REPRO 1      /* 0 = default, others used for comparison      */

#include <iostream>
#include "vegas.h"
using namespace std;

                    /* copy of (*ctl).fcns                          */

double gfsr_rand(unsigned int w[], int *k);
unsigned int kw_rand(void);
void gfsr_init(long seed);

#endif