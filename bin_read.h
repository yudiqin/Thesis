#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <vector>

// typedef float float32;
//typedef complex<double>   data_c;
typedef double            data_dd;
typedef complex<double>            data_cd;

std::vector<data_dd> readFile(const char *filename)
{
    // open the file:
    std::streampos fileSize;
    std::ifstream file(filename, std::ios::binary);

    // get its size:
    file.seekg(0, std::ios::end);
    fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // read the data:
    std::vector<data_dd> fileData(fileSize);
    file.read((char *)&fileData[0], fileSize);
    // variable     type
    // fileData     std::vector<float32>
    // fileData[0]  float32
    // &fileData[0] float32*
    // (char*)&fileData[0]  char*

    return fileData;
}



