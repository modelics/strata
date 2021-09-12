// Copyright 2021 Shashwat Sharma and Piero Triverio

// This file is part of Strata.

// Strata is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Strata is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Strata.  If not, see <https://www.gnu.org/licenses/>.

/************************ lapacke_interface.hpp ************************

 * Wrapper functions for quick access to Lapacke routines.
 *
 * Author: Shashwat Sharma
 * Created on: December 29, 2019

 ***************************************************************/


#ifndef LPCKINT_H
#define LPCKINT_H

#ifndef CBLAS_LAYOUT
#define CBLAS_LAYOUT CBLAS_ORDER
#endif


#include <complex>
#include <vector>
#include <string>
#include <cstring>


namespace strata
{

	// ====================================================
	// Complex numbers
	// ====================================================

	void ComputeLU(std::vector<std::complex<double>> &A, std::vector<int> &ipiv, int m, int n, bool row_major);

	void ComputeInverse(std::vector<std::complex<double>> &A, int n, bool row_major);

	void ComputeInverseFactored(std::vector<std::complex<double>> &A, std::vector<int> &ipiv, int n, bool row_major);

	void SolveLeastSquares(std::vector<std::complex<double>> &A, std::vector<std::complex<double>> &b, int m, int n, int nrhs, bool row_major);
	
	void ComputeSVD(std::vector<std::complex<double>> &A, int m, int n, std::vector<double> &s, std::vector<std::complex<double>> &u, std::vector<std::complex<double>> &vt, bool row_major);

	void ComputeEigenvalues(std::vector<std::complex<double>> &A, int n, std::vector<std::complex<double>> &w, bool row_major);
	
	void MatMatMult_Cblas(std::vector<std::complex<double>> &A, std::vector<std::complex<double>> &B, std::vector<std::complex<double>> &C, int rows_A, int rows_B, std::complex<double> a, std::complex<double> b, bool row_major);


	// ====================================================
	// Real numbers
	// ====================================================

	void ComputeLU(std::vector<double> &A, std::vector<int> &ipiv, int m, int n, bool row_major);

	void ComputeInverse(std::vector<double> &A, int n, bool row_major);

	void MatMatMult_Cblas(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C, int rows_A, int rows_B, double a, double b, bool row_major);

}


#endif

