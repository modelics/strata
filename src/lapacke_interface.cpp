// Author: Shashwat Sharma

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

/************************ lapacke_interface.cpp ************************

 * Wrapper functions for quick access to Lapacke routines.
 *
 * Author: Shashwat Sharma
 * Created on: December 29, 2019

 ***************************************************************/


#include <iostream>
 
#include <cblas.h>
#include <lapacke.h>

#include "lapacke_interface.hpp"


namespace strata
{

	// ====================================================
	// Complex numbers
	// ====================================================

	/*! \brief Function to LU-factorize a general m x n complex<double> matrix. The input matrix is overwritten with the factors.*/
	void ComputeLU(std::vector<std::complex<double>> &A, std::vector<int> &ipiv, int m, int n, bool row_major = true)
	{

		ipiv.resize(m*n);

		int lda, matrix_layout;
		if (row_major)
		{
			matrix_layout = LAPACK_ROW_MAJOR;
			lda = n;
		}
		else
		{
			matrix_layout = LAPACK_COL_MAJOR;
			lda = m;
		}

		int info = LAPACKE_zgetrf(matrix_layout, m, n, (lapack_complex_double *) &A[0], lda, &ipiv[0]);

		return;

	}


	/*! \brief Function to factorize and invert a square complex<double> matrix. The input matrix is overwritten with its inverse.*/
	void ComputeInverse(std::vector<std::complex<double>> &A, int n, bool row_major = true)
	{

		int matrix_layout = LAPACK_ROW_MAJOR;

		if (!row_major)
			matrix_layout = LAPACK_COL_MAJOR;

		std::vector<int> ipiv;
		ComputeLU(A, ipiv, n, n, row_major);

		int lda = n;

		int info = LAPACKE_zgetri(matrix_layout, n, (lapack_complex_double *) &A[0], lda, &ipiv[0]);

		return;

	}


	/*! \brief Function to invert a factored square complex<double> matrix. The input factored matrix is overwritten with its inverse.*/
	void ComputeInverseFactored(std::vector<std::complex<double>> &A, std::vector<int> &ipiv, int n, bool row_major = true)
	{

		int matrix_layout = LAPACK_ROW_MAJOR;

		if (!row_major)
			matrix_layout = LAPACK_COL_MAJOR;
	
		int lda = n;
		int info = LAPACKE_zgetri(matrix_layout, n, (lapack_complex_double *) &A[0], lda, &ipiv[0]);

		return;

	}


	/*! \brief Function to solve a system of equations using minimum-norm least squares via an SVD of the input matrix. Note: RHS vector is overwritten with the solution.*/
	void SolveLeastSquares(std::vector<std::complex<double>> &A, std::vector<std::complex<double>> &b, int m, int n, int nrhs, bool row_major = true)
	{

		int lda, ldb, matrix_layout;
		if (row_major)
		{
			matrix_layout = LAPACK_ROW_MAJOR;
			lda = n;
			ldb = nrhs;
		}
		else
		{
			matrix_layout = LAPACK_COL_MAJOR;
			lda = m;
			ldb = std::max(m, n);
		}

		std::vector<double> s (std::min(m, n));
		double rcond = -1.0;
		int rank;
		
		int info = LAPACKE_zgelss(matrix_layout, m, n, nrhs, (lapack_complex_double *) &A[0], lda, (lapack_complex_double *) &b[0], ldb, &s[0], rcond, &rank);

		if (info < 0)
			std::cout << "[WARNING] SolveLeastSquares(): LAPACKE_zgelss returned " << info << "; parameter " << std::abs(info) << " had an illegal value." << std::endl;
		else if (info > 0)
			std::cout << "[WARNING] SolveLeastSquares(): LAPACKE_zgelss returned " << info << "; did not converge." << std::endl;

		return;

	}


	/*! \brief Function to compute the SVD of a matrix, with the vectors u and vt.*/
	void ComputeSVD(std::vector<std::complex<double>> &A, int m, int n, std::vector<double> &s, std::vector<std::complex<double>> &u, std::vector<std::complex<double>> &vt, bool row_major = true)
	{

		int matrix_layout = LAPACK_ROW_MAJOR;
		int lda = n, ldu = m, ldvt = n;
		u.resize(ldu*m);
		vt.resize(ldvt*n);
		
		if (!row_major)
		{
			lda = m;
			matrix_layout = LAPACK_COL_MAJOR;
		}

		std::vector<double> superb (std::min(m, n));
		s.resize(std::min(m, n));

		int info = LAPACKE_zgesvd(matrix_layout, 'A', 'A', m, n, 
								  (lapack_complex_double *) &A[0], lda, &s[0],
								  (lapack_complex_double *) &u[0], ldu, 
								  (lapack_complex_double *) &vt[0], ldvt, &superb[0]);

		if (info < 0)
			std::cout << "[WARNING] ComputeSVD(): LAPACKE_zgesvd returned " << info << "; parameter " << std::abs(info) << " had an illegal value." << std::endl;
		else if (info > 0)
			std::cout << "[WARNING] ComputeSVD(): LAPACKE_zgesvd returned " << info << "; did not converge." << std::endl;

		return;

	}


	/*! \brief Function to compute the eigenvalues of a square matrix.*/
	void ComputeEigenvalues(std::vector<std::complex<double>> &A, int n, std::vector<std::complex<double>> &w, bool row_major = true)
	{

		int matrix_layout = LAPACK_ROW_MAJOR;

		if (!row_major)
			matrix_layout = LAPACK_COL_MAJOR;

		int lda = n, ldvl = n, ldvr = n;
		std::complex<double> *vl = NULL, *vr = NULL;	
		w.resize(n);

		int info = LAPACKE_zgeev(matrix_layout, 'N', 'N', n, (lapack_complex_double *) &A[0], (lapack_int) lda, (lapack_complex_double *) &w[0], (lapack_complex_double *) vl, ldvl, (lapack_complex_double *) vr, ldvr);

		if (info < 0)
			std::cout << "[WARNING] ComputeEigenvalues(): LAPACKE_zgeev returned " << info << "; parameter " << std::abs(info) << " had an illegal value." << std::endl;
		else if (info > 0)
			std::cout << "[WARNING] ComputeEigenvalues(): LAPACKE_zgeev returned " << info << "; did not converge." << std::endl;

		return;

	}

	
	/*! \brief Function to multiply two matrices of compatible dimensions.*/
	void MatMatMult_Cblas(std::vector<std::complex<double>> &A, std::vector<std::complex<double>> &B, std::vector<std::complex<double>> &C, int rows_A, int rows_B, std::complex<double> a = 1.0, std::complex<double> b = 0.0, bool row_major = true)
	{

		int m = rows_A;
		int n = B.size()/rows_B;
		int k = rows_B;
		int lda = k;
		int ldb = n;
		int ldc = n;

		C.resize(m*n);
		
		CBLAS_LAYOUT matrix_layout = CblasRowMajor;
		
		if (!row_major)
		{
			matrix_layout = CblasColMajor;
			lda = m;
			ldb = k;
			ldc = m;
		}

		cblas_zgemm(matrix_layout, CblasNoTrans, CblasNoTrans, m, n, k, (lapack_complex_double *) &a, (lapack_complex_double *) &A[0], lda, (lapack_complex_double *) &B[0], ldb, (lapack_complex_double *) &b, (lapack_complex_double *) &C[0], ldc);
		// cblas_zgemm(matrix_layout, CblasNoTrans, CblasNoTrans, m, n, k, (double *) &a, (double *) &A[0], lda, (double *) &B[0], ldb, (double *) &b, (double *) &C[0], ldc);


		return;

	}


	// ====================================================
	// Real numbers
	// ====================================================
	
	/*! \brief Function to LU-factorize a general m x n double matrix. The input matrix is overwritten with the factors.*/
	void ComputeLU(std::vector<double> &A, std::vector<int> &ipiv, int m, int n, bool row_major = true)
	{

		ipiv.resize(m*n);

		int lda, matrix_layout;
		if (row_major)
		{
			matrix_layout = LAPACK_ROW_MAJOR;
			lda = n;
		}
		else
		{
			matrix_layout = LAPACK_COL_MAJOR;
			lda = m;
		}

		int info = LAPACKE_dgetrf(matrix_layout, m, n, &A[0], lda, &ipiv[0]);

		return;

	}


	/*! \brief Function to factorize and invert a square double matrix. The input matrix is overwritten with its inverse.*/
	void ComputeInverse(std::vector<double> &A, int n, bool row_major = true)
	{

		int matrix_layout = LAPACK_ROW_MAJOR;

		if (!row_major)
			matrix_layout = LAPACK_COL_MAJOR;

		std::vector<int> ipiv;
		ComputeLU(A, ipiv, n, n, row_major);

		int lda = n;

		int info = LAPACKE_dgetri(matrix_layout, n, &A[0], lda, &ipiv[0]);

		return;

	}

	/*! \brief Function to multiply two matrices of compatible dimensions.*/
	void MatMatMult_Cblas(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C, int rows_A, int rows_B, double a = 1.0, double b = 0.0, bool row_major = true)
	{

		int m = rows_A;
		int n = B.size()/rows_B;
		int k = rows_B;
		int lda = k;
		int ldb = n;
		int ldc = n;

		C.resize(m*n);
		
		CBLAS_LAYOUT matrix_layout = CblasRowMajor;
		
		if (!row_major)
		{
			matrix_layout = CblasColMajor;
			lda = m;
			ldb = k;
			ldc = m;
		}

		cblas_dgemm(matrix_layout, CblasNoTrans, CblasNoTrans, m, n, k, a, &A[0], lda, &B[0], ldb, b, &C[0], ldc);

		return;

	}

}


