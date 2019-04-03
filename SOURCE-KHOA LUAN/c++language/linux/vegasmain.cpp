/**
 * Created by PhpStorm.
 * User: Huynh Thi Ha Vy
 * Date: 03/04/2019
 * Time: 4:24 PM
 */

#include <iostream>
#include <cmath>
using namespace std;
#include "vegas.h"

#define DIMENSION 5
#define FUNCTIONS 4
#define RANDTEST 314159
#ifndef PI
#define PI     3.14159265358979323846
#endif

int t_gfsr_k;
unsigned int t_gfsr_m[SR_P];
double gfsr_norm;

#define WIDTH 0.2

int main(){

    
    gfsr_init(54) ; 
}