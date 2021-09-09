
#ifndef FILEIO_H
#define FILEIO_H

#include <complex>
#include <string.h>
#include <vector>


const int IO_ROW_MAJOR = 1;
const int IO_COL_MAJOR = 2;


void PrintMatToFile_ASCII(std::complex<double> *mat, int rows, int cols, std::string filename, std::complex<double> scale = 1.0, int ordering = IO_ROW_MAJOR);
void PrintVecToFile_ASCII(std::complex<double> *vec, int rows, std::string filename, std::complex<double> scale = 1.0);

void PrintMatToFile_ASCII(std::vector<std::complex<double>> &mat, int rows, int cols, std::string filename, std::complex<double> scale = 1.0, int ordering = IO_ROW_MAJOR);
void PrintVecToFile_ASCII(std::vector<std::complex<double>> &vec, int rows, std::string filename, std::complex<double> scale = 1.0);

void ExportTable_ASCII(std::vector<double *> &table, std::vector<std::string> &headers, int rows, std::string filename);
void ExportTable_ASCII(std::vector<std::complex<double> *> &table, std::vector<std::string> &headers, int rows, std::string filename);


void PrintMatToFile_HDF5(std::complex<double> *mat, int rows, int cols, std::string filename, std::complex<double> scale = 1.0, int ordering = IO_ROW_MAJOR);
void PrintVecToFile_HDF5(std::complex<double> *vec, int rows, std::string filename, std::complex<double> scale = 1.0);

void PrintMatToFile_HDF5(std::vector<std::complex<double>> &mat, int rows, int cols, std::string filename, std::complex<double> scale = 1.0, int ordering = IO_ROW_MAJOR);
void PrintVecToFile_HDF5(std::vector<std::complex<double>> &vec, int rows, std::string filename, std::complex<double> scale = 1.0);


#endif
