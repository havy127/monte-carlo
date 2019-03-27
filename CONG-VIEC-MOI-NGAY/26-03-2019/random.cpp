#include <iostream>
using namespace std;

//preprocess
#define SR_P 1279
#define SR_Q 418
//Globle variable
unsigned int gfsr_m[SR_P];
double gfsr_normn;
double xrand; 
int gfsr_k; 
// declare function 
double gfsr_rand(unsigned int w[], int *k);
unsigned int kw_rand(void);
void gfsr_init(long seed);

int main(){
    
    //xrand =  gfsr_rand(gfsr_m, &gfsr_k);
    //cout << xrand << endl;
    cout << "&gfsr_k : " << &gfsr_k << endl;
    cout << "(*k)++ : " << (&gfsr_k) +1  << endl;
    return 0;
}

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
