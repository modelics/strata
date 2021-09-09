

#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>


#include "fileIO.hpp"


void PrintMatToFile_ASCII(std::complex<double> *mat, int rows, int cols, std::string filename, std::complex<double> scale, int ordering)
{

	std::ofstream file (filename);
	file.precision(15);

	if (ordering == IO_ROW_MAJOR)
	{

		for (int ii = 0; ii < rows; ii++)
		{
			for (int jj = 0; jj < cols; jj++)
			{
				file << std::real(mat[ii*cols + jj]*scale) << ", " << std::imag(mat[ii*cols + jj]*scale) << ", " << std::flush;
			}
			file << std::endl;
		}
		
	}
	else if (ordering == IO_COL_MAJOR)
	{

		for (int jj = 0; jj < cols; jj++)
		{
			for (int ii = 0; ii < rows; ii++)
			{
				file << std::real(mat[ii + rows*jj]*scale) << ", " << std::imag(mat[ii + rows*jj]*scale) << ", " << std::flush;
			}
			file << std::endl;
		}

	}
	else
	{
		std::cout << "[ERROR] PrintMatToFile_ASCII(): Invalid matrix ordering requested." << std::endl;
		file.close();
		return;
	}

	std::cout << "[EXPORT] Wrote < " << filename << " > to disk." << std::endl;
	
	file.close();

	return;
	
}


void PrintMatToFile_ASCII(std::vector<std::complex<double>> &mat, int rows, int cols, std::string filename, std::complex<double> scale, int ordering)
{

	std::ofstream file (filename);
	file.precision(15);

	if (ordering == IO_ROW_MAJOR)
	{

		for (int ii = 0; ii < rows; ii++)
		{
			for (int jj = 0; jj < cols; jj++)
			{
				file << std::real(mat[ii*cols + jj]*scale) << ", " << std::imag(mat[ii*cols + jj]*scale) << ", " << std::flush;
			}
			file << std::endl;
		}
		
	}
	else if (ordering == IO_COL_MAJOR)
	{

		for (int jj = 0; jj < cols; jj++)
		{
			for (int ii = 0; ii < rows; ii++)
			{
				file << std::real(mat[ii + rows*jj]*scale) << ", " << std::imag(mat[ii + rows*jj]*scale) << ", " << std::flush;
			}
			file << std::endl;
		}

	}
	else
	{
		std::cout << "[ERROR] PrintMatToFile_ASCII(): Invalid matrix ordering requested." << std::endl;
		file.close();
		return;
	}

	std::cout << "[EXPORT] Wrote < " << filename << " > to disk." << std::endl;
	
	file.close();

	return;
	
}


void ExportTable_ASCII(std::vector<double *> &table, std::vector<std::string> &headers, int rows, std::string filename)
{

	std::ofstream file (filename);
	file.precision(15);

	// Print headers
	file << headers[0] << std::flush;
	for (int ii = 1; ii < headers.size(); ii++)
		file << "; " << headers[ii] << std::flush;
	file << std::endl;

	// Print table
	for (int jj = 0; jj < rows; jj++)
	{
		file << table[0][jj] << std::flush;
		for (int ii = 1; ii < table.size(); ii++)
		{
			file << " " << table[ii][jj] << std::flush;
		}
		file << std::endl;
	}

	std::cout << "[EXPORT] Wrote < " << filename << " > to disk." << std::endl;

	return;

}


void ExportTable_ASCII(std::vector<std::complex<double> *> &table, std::vector<std::string> &headers, int rows, std::string filename)
{

	std::ofstream file (filename);
	file.precision(15);

	// Print headers
	file << headers[0] << std::flush;
	for (int ii = 1; ii < headers.size(); ii++)
		file << "; " << headers[ii] << std::flush;
	file << std::endl;

	// Print table
	for (int jj = 0; jj < rows; jj++)
	{
		file << table[0][jj] << std::flush;
		for (int ii = 1; ii < table.size(); ii++)
		{
			file << " " << std::real(table[ii][jj]) << " " << std::imag(table[ii][jj]) << std::flush;
		}
		file << std::endl;
	}

	std::cout << "[EXPORT] Wrote < " << filename << " > to disk." << std::endl;

	return;

}
