#include "BVector_Filter_OneSample_Sub.h"


void Top(data_c *tb_in, data_c *tb_out, data_d *coeff, int N)
{
	const int NUM_ROW_ = NUM_ROW;
	const int NUM_COLUM_ = NUM_COLUM;
	const int NUM_COEFF_ = NUM_COEFF;
	const int depth_factor = 1;

#pragma HLS INTERFACE mode=ap_memory port=tb_in depth=(NUM_ROW_*NUM_COLUM_)*depth_factor
#pragma HLS INTERFACE mode=ap_memory port=tb_out depth=((NUM_ROW_*4)*NUM_COLUM_)*depth_factor
#pragma HLS INTERFACE mode=ap_memory port=coeff depth=NUM_COEFF_
//#pragma HLS INTERFACE mode=m_axi port=tb_in depth=(NUM_ROW_*NUM_COLUM_)*depth_factor bundle=tb_in_bundle
//#pragma HLS INTERFACE mode=m_axi port=tb_out depth=((NUM_ROW_*4)*NUM_COLUM_)*depth_factor bundle=tb_out_bundle
//#pragma HLS INTERFACE mode=m_axi port=coeff depth=NUM_COEFF_ bundle=coeff_bundle


/* init B*/
	data_d B[NUM_COEFF];
#pragma HLS ARRAY_PARTITION variable=B type=complete dim=1
read_coeff_loop:
    for (int i = 0; i < NUM_COEFF; i++)
    {
#pragma HLS pipeline II = 1 rewind
        B[i] = coeff[i];
    }
/* end init B*/


/* main */
//	loop_compute((data_c *)tb_in, B, N, (data_c *)tb_out);
#pragma HLS DATAFLOW
loop_compute_row:
    for (int j = 0; j < NUM_COLUM; j++)
    {
	loop_compute_column:
    	for (int i = 0; i < N*4; i++)
        {
            data_c Y;
            data_c input_temp = tb_in[j + i / 4 * NUM_COLUM];

//
 			data_d input_real = input_temp.real() * 2.0;
            data_d input_imag = input_temp.imag() * 2.0;
#pragma HLS LOOP_TRIPCOUNT  min=61440*4 max=61440*4 avg=61440*4
//#pragma HLS DATAFLOW
            Fir_real(input_real, B, &Y, i);
            Fir_imag(input_imag, B, &Y, i);
            writeBackY(Y, i, j, (data_c*)tb_out);
        }
    }
/* end main */

}


void Fir_real(data_d input_real, data_d B[NUM_COEFF], data_c *Y, int i)
{
	static int oldest_Z_idx_real = 0;
    const int Z_idx_real = oldest_Z_idx_real;
    static data_d Z_real[NUM_COEFF-1];  // NUM_COEFF = 257
#pragma HLS ARRAY_PARTITION variable=Z_real type=complete

    data_d Z_buffer_real[NUM_COEFF-1];
#pragma HLS ARRAY_PARTITION variable=Z_buffer_real type=complete

    data_d X_real;
    data_d zeros_real = 0.0;

#pragma HLS PIPELINE II=1
//#pragma HLS DATAFLOW

//used to clear the Z and index when start from a new column
    if (i==0)
    {
		for(int i = 0; i<NUM_COEFF - 1; i++)
		{
#pragma HLS UNROLL  // max=256
			Z_real[i] = 0;
		}
		oldest_Z_idx_real = 0;
    }

//upsample
    if(i % 4 == 0)
        X_real = input_real;
    else
        X_real = zeros_real;

    copyUpdateZ_real(X_real, Z_real, oldest_Z_idx_real, Z_idx_real, Z_buffer_real);

    // y = f1(b, x)
    vectorInit_real(B, X_real, Y);

    // y = f2(b, z, y)
    multiAccumulateCalc_real(B, Z_idx_real, Y, Z_buffer_real);
}

void Fir_imag(data_d input_imag, data_d B[NUM_COEFF], data_c *Y, int i)
{
	static int oldest_Z_idx_imag = 0;
    const int Z_idx_imag = oldest_Z_idx_imag;
    static data_d Z_imag[NUM_COEFF-1];  // NUM_COEFF = 257
#pragma HLS ARRAY_PARTITION variable=Z_imag type=complete

    data_d Z_buffer_imag[NUM_COEFF-1];
#pragma HLS ARRAY_PARTITION variable=Z_buffer_imag type=complete

    data_d X_imag;
    data_d zeros_imag = 0.0;

#pragma HLS PIPELINE II=1
//#pragma HLS DATAFLOW
//used to clear the Z and index when start from a new column
    if (i==0)
    {
		for(int i = 0; i<NUM_COEFF - 1; i++)
		{
#pragma HLS UNROLL  // max=256
			Z_imag[i] = 0;
		}
		oldest_Z_idx_imag = 0;
    }
//upsample here
    if(i % 4 == 0)
        X_imag = input_imag;
    else
        X_imag = zeros_imag;

    copyUpdateZ_imag(X_imag, Z_imag, oldest_Z_idx_imag, Z_idx_imag, Z_buffer_imag);

    // y = f1(b, x)
    vectorInit_imag(B, X_imag, Y);

    // y = f2(b, z, y)
    multiAccumulateCalc_imag(B, Z_idx_imag, Y, Z_buffer_imag);
}


void copyUpdateZ_real(data_d X_real, data_d Z_real[NUM_COEFF - 1], int &oldest_Z_idx_real, const int &Z_idx_real, data_d Z_buffer_real[NUM_COEFF - 1])
{
#pragma HLS INLINE
	for(int i = 0; i < NUM_COEFF - 1; i++)
	{
		Z_buffer_real[i] = Z_real[i];
	}
	Z_real[Z_idx_real] = X_real;
	oldest_Z_idx_real = (Z_idx_real + 1) % (NUM_COEFF - 1);
}


void copyUpdateZ_imag(data_d X_imag, data_d Z_imag[NUM_COEFF - 1], int &oldest_Z_idx_imag, const int &Z_idx_imag, data_d Z_buffer_imag[NUM_COEFF - 1])
{
#pragma HLS INLINE
	for(int i = 0; i < NUM_COEFF - 1; i++)
	{
		Z_buffer_imag[i] = Z_imag[i];
	}
	Z_imag[Z_idx_imag] = X_imag;
	oldest_Z_idx_imag = (Z_idx_imag + 1) % (NUM_COEFF - 1);
}

void vectorInit_real(data_d B[NUM_COEFF], data_d X_real, data_c *Y)
{
#pragma HLS INLINE
	VEC_INIT:
    // Y[j] = X[j] * B[0];
    Y->real(X_real * B[0]);
}


void multiAccumulateCalc_real(data_d B[NUM_COEFF], const int &Z_idx_real, data_c *Y, data_d Z_buffer_real[NUM_COEFF - 1])
{
#pragma HLS INLINE
    const int unroll_factor = 4;
    data_d part_sum_real;
    MAC:
    for (int i = NUM_COEFF - 1; i != 0; i--)
    {
   	    const int current_Z_idx_real = (Z_idx_real- i + (NUM_COEFF-1)) % (NUM_COEFF-1);
        part_sum_real += Z_buffer_real[current_Z_idx_real] * B[i];
		if(i % unroll_factor == 1)  // i= 4,3,2,1 -> i%4= 0,3,2,1; i%2= 0,1,0,1
		{
             Y->real(Y->real() + part_sum_real);
             part_sum_real = 0;
		}
    }
}


void vectorInit_imag(data_d B[NUM_COEFF], data_d X_imag, data_c *Y)
{
#pragma HLS INLINE
	VEC_INIT:
    // Y[j] = X[j] * B[0];
    Y->imag(X_imag * B[0]);

}

void multiAccumulateCalc_imag(data_d B[NUM_COEFF], const int &Z_idx_imag, data_c *Y, data_d Z_buffer_imag[NUM_COEFF - 1])
{
#pragma HLS INLINE
    const int unroll_factor = 8;
    data_d part_sum_imag;
    MAC:
    for (int i = NUM_COEFF - 1; i != 0; i--)
    {
    	const int current_Z_idx_imag = (Z_idx_imag- i + (NUM_COEFF-1)) % (NUM_COEFF-1);
        part_sum_imag += Z_buffer_imag[current_Z_idx_imag] * B[i];
		if(i % unroll_factor == 1)  // i= 4,3,2,1 -> i%4= 0,3,2,1; i%2= 0,1,0,1
		{
             Y->imag(Y->imag() + part_sum_imag);
             part_sum_imag = 0;
		}
    }
}


void writeBackY(data_c Y, int i, int j, data_c *tb_out)
{
//#pragma HLS INLINE
	tb_out[i + j * NUM_ROW * 4] = Y;
}



