#include "BVector_Filter_SubFir.h"

// extern "C" {

void Top(data_c *tb_in, data_c *tb_out, data_d *coeff, const int &N)
{


#pragma HLS INTERFACE axis port = tb_in
#pragma HLS INTERFACE axis port = tb_out
#pragma HLS INTERFACE axis port = coeff

    data_c zeros[NUM_COLUM]; // NUM_COLUM==32

    // init: read the coefficient as local variable
    data_d B[NUM_COEFF];
#pragma HLS ARRAY_PARTITION variable=B type=complete dim=1
#pragma HLS STREAM variable=zeros depth=3
#pragma HLS STREAM variable=B depth=3

//    data_c X[NUM_COLUM]; // NUM_COLUM==32
//    data_c Y[NUM_COLUM];
//#pragma HLS ARRAY_PARTITION variable = X type = complete
//#pragma HLS ARRAY_PARTITION variable = Y type = complete

//    static data_c Z[NUM_COLUM][NUM_COEFF - 1]; // [32][256]
      data_c Z[NUM_COLUM][NUM_COEFF - 1]; // [32][256]
#pragma HLS ARRAY_PARTITION variable=Z type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Z type=cyclic factor=4 dim=2
//    static int oldest_Z_index = 0;
       int oldest_Z_index = 0;

#pragma HLS DATAFLOW

    read_coff:
    for (int i = 0; i < NUM_COEFF; i++)
    {
#pragma HLS pipeline II = 1 rewind
        B[i] = coeff[i];
    }

    loop_compute_column:
    for (int i = 0; i < N; i++) // N==61440 is the number of row. so the operation focus on every row
    {
        data_c X[NUM_COLUM]; // NUM_COLUM==32
        data_c Y[NUM_COLUM];
#pragma HLS ARRAY_PARTITION variable=X type=complete
#pragma HLS ARRAY_PARTITION variable=Y type=complete

#pragma HLS LOOP_TRIPCOUNT min=61440 max=61440 avg=61440
#pragma HLS DATAFLOW

        readInX(X, (data_c *)tb_in, i);

        // 2.Fir(0)
        Fir(X, Y, B, Z, oldest_Z_index);
//        UpdateZindex(oldest_Z_index);
        writeBackY(Y, (data_c *)tb_out, i, 0);

        //Fir(1)
        Fir(zeros, Y, B, Z, oldest_Z_index);
//        UpdateZindex(oldest_Z_index);
        writeBackY(Y, (data_c *)tb_out, i, 1);

        //Fir(2)
        Fir(zeros, Y, B, Z, oldest_Z_index);
//        UpdateZindex(oldest_Z_index);
        writeBackY(Y, (data_c *)tb_out, i, 2);

        //Fir(3)
        Fir(zeros, Y, B, Z, oldest_Z_index);
//        UpdateZindex(oldest_Z_index);
        writeBackY(Y, (data_c *)tb_out, i, 3);


    }
}

// 1.read in
void readInX(data_c X[NUM_COLUM], data_c *tb_in, int i)
{
    for (int j = 0; j < NUM_COLUM; j++)
    {
#pragma HLS PIPELINE II = 1 rewind
        // X= 2*tb_in
        X[j] = tb_in[j + i * NUM_COLUM] * 2.0; // take all elements from each row of the tb_in
    }
}

 //Compute Fir function
void writeBackY(data_c Y[NUM_COLUM], data_c *tb_out, int i, int k)
{
    for (int j = 0; j < NUM_COLUM; j++)
    {
#pragma HLS pipeline II = 1 rewind
        tb_out[i * 4 + k + j * NUM_ROW * 4] = Y[j];
    }
}

void Fir(data_c X[NUM_COLUM], data_c Y[NUM_COLUM], data_d B[NUM_COEFF], data_c Z[NUM_COLUM][NUM_COEFF - 1], int &oldest_Z_index) // NUM_COLUM==32, NUM_COEFF==257
{
    // const int unroll_factor = 16;

#pragma HLS ARRAY_PARTITION variable = X type = complete dim = 1
#pragma HLS ARRAY_PARTITION variable = Y type = complete dim = 1
#pragma HLS ARRAY_PARTITION variable = B type = complete dim = 1

   const int Z_index = oldest_Z_index;

//    static data_c Z[NUM_COLUM][NUM_COEFF - 1]; // [32][256]
//#pragma HLS ARRAY_PARTITION variable = Z type = complete dim = 1
//#pragma HLS ARRAY_PARTITION variable = Z type = cyclic factor = 4 dim = 2
//    static int oldest_Z_index = 0;

    // data_c part_sum[NUM_COLUM];
    // #pragma HLS ARRAY_PARTITION variable = part_sum type = complete

#pragma HLS DATAFLOW
    // y = f1(b, x)
    vectorInit(B, X, Y);

    // y = f2(b, z, y)
    multiAccumulateCalc(B, Z, Z_index, Y);

    // z = f3(x, z, oldest_Z_index)
//    updateZ(X, Z, oldest_Z_index, Z_index);
	WriteToZ(X, Z, Z_index);
	UpdateZindex(Z_index,oldest_Z_index);

}

/*

*/

// y = f1(b, x)
void vectorInit(data_d B[NUM_COEFF], data_c X[NUM_COLUM], data_c Y[NUM_COLUM])
{
VEC_INIT:
    for (int j = 0; j < NUM_COLUM; j++)
    {
#pragma HLS pipeline II = 1 rewind
//#pragma HLS UNROLL
        Y[j].real(X[j].real() * B[0]);
        Y[j].imag(X[j].imag() * B[0]);
    }
}

// y = f2(b, z, y)
void multiAccumulateCalc(data_d B[NUM_COEFF], data_c Z[NUM_COLUM][NUM_COEFF - 1], const int &Z_index, data_c Y[NUM_COLUM])
{
#pragma HLS DATAFLOW
    const int unroll_factor = 8;
    data_c part_sum[NUM_COLUM];
#pragma HLS ARRAY_PARTITION variable = part_sum type = complete


MAC:
//    for (int i = NUM_COEFF - 1; i != 0; i--)
	  for (int i = 1; i < NUM_COEFF; i++)
    {
#pragma HLS pipeline II = 4 rewind// single complex=4, two complex=2
#pragma HLS unroll factor = unroll_factor

        // int i_complement = (NUM_COEFF - 1) - i;
        const int current_Z_index = (Z_index + (NUM_COEFF - 1) - i) % (NUM_COEFF - 1);

    VEC_MAC:
        for (int j = 0; j < NUM_COLUM; j++)
        {
            // part_sum[j] += Z[j][current_Z_index] * B[i];
            part_sum[j].real(part_sum[j].real() + Z[j][current_Z_index].real() * B[i]);
            part_sum[j].imag(part_sum[j].imag() + Z[j][current_Z_index].imag() * B[i]);

            if (i % unroll_factor == 1) // i= 4,3,2,1 -> i%4= 0,3,2,1; i%2= 0,1,0,1
            {
                Y[j].real(Y[j].real() + part_sum[j].real());
                Y[j].imag(Y[j].imag() + part_sum[j].imag());

                part_sum[j].real(0.0);
                part_sum[j].imag(0.0);
            }
        }
    }
}

// z = f3(x, z, oldest_Z_index)
//void updateZ(data_c X[NUM_COLUM], data_c Z[NUM_COLUM][NUM_COEFF - 1], int &oldest_Z_index, const int &Z_index)
//{
//#pragma HLS DATAFLOW
////Write_to_Z:
////    for (int j = 0; j < NUM_COLUM; j++)
////    {
////#pragma HLS pipeline II = 1 rewind
//////#pragma HLS unroll
////        Z[j][Z_index] = X[j];
////    }
////    oldest_Z_index = (oldest_Z_index + 1) % (NUM_COEFF - 1);
//	WriteToZ(X, Z, Z_index);
//	UpdateZindex(Z_index,oldest_Z_index);
//
//}

void WriteToZ(data_c X[NUM_COLUM], data_c Z[NUM_COLUM][NUM_COEFF - 1], const int &Z_index)
{
		Write_to_Z:
	    for (int j = 0; j < NUM_COLUM; j++)
	    {
	#pragma HLS pipeline II = 1 rewind
	//#pragma HLS unroll
	        Z[j][Z_index] = X[j];
	    }

}

void UpdateZindex(const int &Z_index,int &oldest_Z_index)
{

	oldest_Z_index = (Z_index + 1) % (NUM_COEFF - 1);

}
