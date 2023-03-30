#ifndef __BVECTOR_FILTER_H__
#define __BVECTOR_FILTER_H__

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include "math.h"
#include <cstdlib>
#include <string>

using namespace std;

//the input matrix is the initial matrix has been tranposed
#define NUM_ROW  61440 //6
#define NUM_COLUM  32 //3
#define NUM_COEFF  257 //5
#define ERR 0.000001

//to define the number of unroll and cyclic factor
#ifndef MAC_unroll_factor
#define MAC_unroll_factor 1
#endif

#ifndef Z_cyclic_factor
#define Z_cyclic_factor 1
#endif

#ifndef SHIFT_unroll_factor
#define SHIFT_unroll_factor 1
#endif

typedef complex<double>   data_c;
//typedef double           data_c;
typedef double            data_d;

//extern "C" {
// Top function
void Top(data_c *tb_in, data_c *tb_out, data_d *coeff, const int& N);
//}

//Sub Function
void readInX(data_c X[NUM_COLUM], data_c *tb_in, int i);
void Fir(data_c X[NUM_COLUM], data_c Y[NUM_COLUM], data_d B[NUM_COEFF], data_c Z[NUM_COLUM][NUM_COEFF - 1], int &oldest_Z_index);
//void UpdateZindex(int &Z_index,int &oldest_Z_index);
void writeBackY(data_c Y[NUM_COLUM], data_c *tb_out, int i, int k);

//sub of Fir
void vectorInit(data_d B[NUM_COEFF], data_c X[NUM_COLUM], data_c Y[NUM_COLUM]);
void multiAccumulateCalc(data_d B[NUM_COEFF], data_c Z[NUM_COLUM][NUM_COEFF - 1], const int &Z_index, data_c Y[NUM_COLUM]);
//void updateZ(data_c X[NUM_COLUM], data_c Z[NUM_COLUM][NUM_COEFF - 1], int &oldest_Z_index,const int &Z_index);
void WriteToZ(data_c X[NUM_COLUM], data_c Z[NUM_COLUM][NUM_COEFF - 1], const int &Z_index);
void UpdateZindex(const int &Z_index, int &oldest_Z_index);


 #endif
