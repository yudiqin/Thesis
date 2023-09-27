#ifndef __BVECTOR_FILTER_H__
#define __BVECTOR_FILTER_H__

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <math.h>
#include <cstdlib>
#include <string>

#include <stdlib.h>
#include <string.h>


using namespace std;

//the input matrix is the initial matrix has been tranposed
#define NUM_ROW  61440
#define NUM_COLUM  32
#define NUM_COEFF  257
#define ERR 0.000001
#define ERR_2 0.00000001


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

//typedef complex<double>   data_c;
//typedef double           data_d;

typedef complex<float>   data_c;
typedef float           data_d;


// Top function
void Top(data_c *tb_in, data_c *tb_out, data_d *coeff, int N);

//Sub Function


//Sub loop_c_c
void Fir_real(data_d input_real, data_d B[NUM_COEFF], data_c *Y, int i);
void Fir_imag(data_d input_imag, data_d B[NUM_COEFF], data_c *Y, int i);
void writeBackY(data_c Y, int i, int j, data_c *tb_out);

//
void copyUpdateZ_real(data_d X_real, data_d Z_real[NUM_COEFF - 1], int &oldest_Z_idx_real, const int &Z_idx_real, data_d Z_buffer_real[NUM_COEFF - 1]);
void copyUpdateZ_imag(data_d X_imag, data_d Z_imag[NUM_COEFF - 1], int &oldest_Z_idx_imag, const int &Z_idx_imag, data_d Z_buffer_imag[NUM_COEFF - 1]);
void vectorInit_real(data_d B[NUM_COEFF], data_d X_real, data_c *Y);
void multiAccumulateCalc_real(data_d B[NUM_COEFF], const int &Z_idx_real, data_c *Y, data_d Z_buffer_real[NUM_COEFF - 1]);
void vectorInit_imag(data_d B[NUM_COEFF], data_d X_imag, data_c *Y);
void multiAccumulateCalc_imag(data_d B[NUM_COEFF], const int &Z_idx_imag, data_c *Y, data_d Z_buffer_imag[NUM_COEFF - 1]);



 #endif
