#include "BVector_Filter_SubFir.h"
#include "bin_read.h"

data_c in_data[NUM_ROW][NUM_COLUM];
data_c out_data[NUM_COLUM][NUM_ROW*4];
data_d in_coeff[NUM_COEFF];
data_c mat_out[NUM_COLUM][NUM_ROW*4];
const int Num = NUM_ROW;

//used to store data from matlab file to compare
data_d data_in_real[NUM_ROW][NUM_COLUM];
data_d data_in_imag[NUM_ROW][NUM_COLUM];
data_d data_out_real[NUM_COLUM][NUM_ROW*4];
data_d data_out_imag[NUM_COLUM][NUM_ROW*4];
data_d in_coeff_temp[1][NUM_COEFF];


int main()
{
     //To gengrate the coeff
//	data_d count = 1;
//	for (int i = 0; i < NUM_COEFF; i++)
//	{
//
//		in_coeff[i] = count;
//		count ++;
//	}
	//To read the coeff
	std::vector<data_d>  coeff_m = readFile("coeff.bin");

	for (int i = 0; i < NUM_COEFF; i++)
	{
			in_coeff[i] = coeff_m[i];
	}

	 // display the coeff
//	cout << "the coeff is" << "\n";
//
//		for(int j=0; j<NUM_COEFF; j++)
//		{
//			std::stringstream ss_r;
//			ss_r << std::fixed << in_coeff[j];
//			cout << ss_r.str() << "\t";
//
//	}
//	cout << endl;

    //read input from file which is same as Matlab input
    std::vector<data_d>  inputTR = readFile("input_real.bin");
    std::vector<data_d>  inputTI = readFile("input_imag.bin");

    //Transpose the bin matrix
    for (int j = 0; j < NUM_COLUM; j++)
	{
		for (int i = 0; i < NUM_ROW; i++)
		{
            
            data_in_real[i][j] = inputTR[i + j * NUM_ROW];
            data_in_imag[i][j] = inputTI[i + j * NUM_ROW];

		}
    }
    //Merge the real and imag part
    for (int i = 0; i < NUM_ROW; i++)
	{
			for (int j = 0; j < NUM_COLUM; j++)
			{

				in_data[i][j].real(data_in_real[i][j]);// = RY[i][j];
				in_data[i][j].imag(data_in_imag[i][j]);//= IY[i][j];
			
			 }
	 }
 // display the input
//         	cout << "the input is"
//         		 << "\n";
//         	for (int i = 0; i < NUM_ROW; i++)
//         	{
//         		for (int j = 0; j < NUM_COLUM; j++)
//         		{
//         			std::stringstream ss_r;
//         			ss_r << std::fixed << in_data[i][j];
//         			cout << ss_r.str() << "\t";
//
//         		}
//         		cout << endl;
//         	}
//         	cout << "\n";

    Top((data_c*)in_data, (data_c*)out_data, in_coeff, Num);//(data_c*)is used to force the matrix to 1D array?

       // display the output
//     	cout << "the output is"
//     		 << "\n";
//     	for (int i = 0; i < NUM_COLUM; i++)
//     	{
//
//     		for (int j = 0; j < NUM_ROW*4; j++)
//     		{
//     			std::stringstream ss_r;
//     			ss_r << std::fixed << out_data[i][j];
//     			cout << ss_r.str() << "\t";
//     		}
//     		cout << endl;
//     	}
//     	cout << "\n";
 

    //read output from  Matlab input
    std::vector<data_d>  outputTR = readFile("output_real.bin");
    std::vector<data_d>  outputTI = readFile("output_imag.bin");

    //Transpose the bin matrix
    for (int j = 0; j < NUM_ROW*4; j++)
	{
		for (int i = 0; i < NUM_COLUM; i++)
		{
            
            data_out_real[i][j] = outputTR[i + j * NUM_COLUM];
            data_out_imag[i][j] = outputTI[i + j * NUM_COLUM];

		}
    }

    //Merge the real and imag part
    for (int i = 0; i < NUM_COLUM ; i++)
	{
		for (int j = 0; j < NUM_ROW*4; j++)
		{

			mat_out[i][j].real(data_out_real[i][j]);
			mat_out[i][j].imag(data_out_imag[i][j]);

		}
	}

    //display the Matlab output
//     	cout << "the Matlab_output is"
//     		 << "\n";
//     	for (int i = 0; i < NUM_COLUM; i++)
//     	{
//     		for (int j = 0; j < NUM_ROW*4; j++)
//     		{
//     			std::stringstream ss_r;
//     			ss_r << std::fixed << mat_out[i][j];
//     			cout << ss_r.str() << "\t";
//
//     		}
//     		cout << endl;
//     	}
//     	cout << "\n";
    //display the Matlab output
//     	cout << "the Matlab_output is"
//     		 << "\n";
// 		for (int j = 0; j < NUM_ROW*4; j++) //Trans_col
//     	{
// 			for (int i = 0; i < NUM_COLUM; i++) //Tranps_row
//     		{
//     			std::stringstream ss_r1;
//     			ss_r1 << std::fixed << mat_out[i][j];
//     			cout << ss_r1.str() << "\t";
//     			std::stringstream ss_r;
//     			ss_r << std::fixed << out_data[i][j];
//     			cout << ss_r.str() << "\t";
//     			printf("row is %d \t column:%d \t",i,j);
//     			cout << "\n";
//
//     		}
//     		cout << endl;
//     	}
//     	cout << "\n";

    //to compare the output and calculatthe error percentage
	int match = 0;
	int mismatch = 0;
	data_d real_error;
	data_d imag_error;

	for (int i = 0; i < NUM_COLUM; i++)
	{
		for (int j = 0; j < NUM_ROW*4; j++)
		{
			//for real part
			if(mat_out[i][j].real() == 0)
				real_error = fabs(out_data[i][j].real());
			else
				real_error = fabs((out_data[i][j].real() - mat_out[i][j].real()) / mat_out[i][j].real());

			//for imag part
			if(mat_out[i][j].imag() == 0)
				imag_error = fabs(out_data[i][j].imag());
			else
			 imag_error = fabs((out_data[i][j].imag() - mat_out[i][j].imag()) / mat_out[i][j].imag());
				

            if((real_error < ERR) && (imag_error < ERR))
				match++;

			else
				mismatch++;
			
		}
		
	}

	printf("OUTPUT match:%d \t mismatch:%d \t",match,mismatch);

	return 0;

}
