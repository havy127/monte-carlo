#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
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

int main(){
    
    //xrand =  gfsr_rand(gfsr_m, &gfsr_k);
    //cout << xrand << endl;
    printf("(*k)++ : %ld /n", (&gfsr_k)++ );
    //cout << "&gfsr_k : " << &gfsr_k << endl;
    //cout << "(*k)++ : " << (*gfsr_k)++ << endl;
    return 0;
}

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