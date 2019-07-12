#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
double random_number_uniform (double a, double b){

   return a + (b - a)*rand()/RAND_MAX;

 }
int main(int argc, char ** argv) 
{ 	int i;
   srand(13); 
   int numSims = 1000000; 
   int numBins = 100; 
   int bins[numBins]; 
   memset(bins, 0x0, sizeof(int) * numBins); // set all the bins to 0 
   const float dist = 10; // 10 km 
   for (i = 0; i < numSims; ++i) { 
   	   float num_random = random_number_uniform(0.0,1.0);
//   	   printf("ok: %f",num_random);
       float diff = (2 * num_random -1) * 5; // random var between -5 and 5 
       int whichBin = (int)(numBins * (diff / 5 + 1) * 0.5); 
       bins[whichBin]++; 
       float time = 30 + diff; 
       float speed = 60 * dist / time; 
   } 
   float sum = 0; 
   for (i = 0; i < numBins; ++i) { 
       float r = bins[i] / (float)numSims; 
       printf("%f %f\n", 5 * (2 * (i /(float)(numBins)) -1), r); 
       sum += r; 
   } 
   fprintf(stderr, "sum %f\n", sum); 
   return 0; 
} 
