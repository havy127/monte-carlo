#include "vegasclass.h"
#include <cmath>
#include <fstream>


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
  w[*k] = w[*k] ^ w[j];
  return((double)w[*k] * gfsr_norm);
}

void randtest()
{
  int i;
  unsigned int rdum = 42;
  double x, x0, sigma; 
  ofstream t_gfsr_m_file;
  t_gfsr_m_file.open ("gfsr_m_file.txt");
  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));

  cout << "(pow(2.0,(long double)(8*sizeof(int))) - 1) : " << gfsr_norm << endl;
  for (i=0; i<SR_P; i++) { 
    rdum = (1812433253*rdum + 314159265);
    t_gfsr_m[i] = rdum;
    t_gfsr_m_file << t_gfsr_m[i] << endl;
    cout << t_gfsr_m[i] << endl;
  }
  t_gfsr_k = -1;
  t_gfsr_m_file.close();
  cout << "Primitive test of the RNG:" << endl;
  cout << " expecting:  mean=0.500000   sigma=0.083333" << endl ; 
  sigma = x0 = 0;
  /*
  * - Tạo ra mãng t_gfsr_m[] chứa các gia trị ngẫu nhiên x (không lớn hơn 1)
  * - x0: là số TB.
  * - Vòng for thứ 2: tính sai số sigma.  
  */
  for (i=0; i< RANDTEST; i++) { 
    x = myrand(t_gfsr_m,&t_gfsr_k);
    // cout << "x=myrand: " << x <<endl ;
    if (fabs(0.5-x)>0.5) cout << "Oh my god, RN=" << x;
    // cout << "x0t=myrand: " << x0 <<endl ;
    x0 += x;
    // cout << "x0=myrand: " << x0 <<endl ;
  }
  cout << "x0=myrand: " << x0 <<endl ;
  x0 = x0 /RANDTEST;
  cout << " calculated: mean=" << x0 ; 

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




