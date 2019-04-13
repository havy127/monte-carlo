#include"vegas.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;

/**================================================== *
 * ==========  variable for RANDTEST  ========== *
 * ================================================== */

#define RANDTEST 314159
int t_gfsr_k;
unsigned int t_gfsr_m[SR_P];
double gfsr_norm;

double myrand(unsigned int w[], int *k)
{
  int j;
  // cout << "myrand with K value: " << *k << endl;
  (*k)++;
  if (*k >= SR_P) *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P) j -= SR_P;
  // cout << "w[*k]: " << w[*k] << endl;
  // cout << "w[j] " << w[j] << endl;
  w[*k] = w[*k] ^ w[j];
  // cout << "w[*k] ^ w[j]: " << w[*k] << endl;
  // cout << "(double)w[*k] * gfsr_norm: " << (double)w[*k] * gfsr_norm << endl;
  return((double)w[*k] * gfsr_norm);
}

void randtest()
{
  int i;
  long seed = time(NULL);
  cout << "time(NULL): " << seed << endl;
  unsigned int rdum = seed;
  double x, x0, sigma; 
  ofstream t_gfsr_m_file, x_numbers_file;
  t_gfsr_m_file.open ("t_gfsr_m_file.txt");
  x_numbers_file.open ("x_numbers_file.txt");

  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));
  cout << "gfsr_norm : " << gfsr_norm << endl;
  cout << "t_gfsr_k : " << &t_gfsr_k << endl;
  for (i=0; i<SR_P; i++) { 
    rdum = (1812433253*rdum + 314159265);
    t_gfsr_m[i] = rdum;
    t_gfsr_m_file << t_gfsr_m[i] << endl;
  }
  t_gfsr_k = -1;
  t_gfsr_m_file.close();
  cout << "Primitive test of the RNG:" << endl;
  cout << " expecting:  mean=0.500000   sigma=0.083333" << endl ; 
  sigma = x0 = 0;
  /*
  * - Tạo ra mãng t_gfsr_m[] chứa các gia trị ngẫu nhiên x (không lớn hơn 1)
  * - x0: là số TB.
  * - Vòng for thứ 2: tính sai số sigma.  RANDTEST
  */ 
  for (i=0; i< RANDTEST; i++) { 

    x = myrand(t_gfsr_m,&t_gfsr_k);
    x_numbers_file << x <<endl;
    add_t_gfsr_k_file << &t_gfsr_k <<endl;
  
    if (fabs(0.5-x)>0.5) cout << "Oh my god, RN=" << x;
    x0 += x;
  }
  x_numbers_file.close();

  x0 = x0 /RANDTEST;

  for (i=0; i < RANDTEST; i++) {
    x = myrand(t_gfsr_m,&t_gfsr_k);
    if (fabs(0.5-x)>0.5) cout << "Oh my god, RN=" << x;
    sigma += (x-x0)*(x-x0); // sai số 
  }
  sigma = sigma / RANDTEST;
  cout << " sigma=" << sigma << "  normalized with " << gfsr_norm ; 
  return;
}

/* =======  End of Section comment block  ======= */



int main (){

    // gfsr_init(12);
    // cout << "gfsr_rand(): " << gfsr_rand() << endl;
    randtest();

    return 0;
}




