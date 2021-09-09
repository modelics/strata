/************************ sommerfeld_integrals.hpp ************************

 * Toolkit for computing Sommerfeld integrals in the complex krho plane.
 * Includes straightforward numerical integration, as well as the
 * Levin-Sidi approach with partition-extrapolation. This is a header-only
 * library.

 * References: 
 * Michalski, Mosig, "Efficient Computation of Sommerfeld Integral Tails -
 * Methods and Algorithms", Jnl. EM Waves & Appl., 2016.
 *
 * Author: Shashwat Sharma
 * Created on: Mar 14, 2020

 *************************************************************************/


#ifndef SI_H
#define SI_H


#include <complex>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp> 
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "libAmosBessel.h"

#include "spectral_MGF.hpp"
#include "bessel_zeros.hpp"


// ==================================================================================
// Declarations
// ==================================================================================

// ====== Computational drivers ======

std::complex<double> SommerfeldIntegrand(SpectralMGF &smgf, double rho, std::complex<double> krho, int component, int order, bool curl = false);
std::complex<double> IntegrateSpectralFarField(SpectralMGF &smgf, double rho, int component, int order, double a, bool curl = false);
std::complex<double> IntegrateSpectralNearField(SpectralMGF &smgf, double rho, int component, int order, double a, bool curl = false);


// ====== Partition-extrapolation integration routines ======

template<class F>
void PartExtrap(const F f, double a, double q, double tol, std::complex<double> &val);

std::complex<double> LevinSidi(int k, std::complex<double> Sk, std::complex<double> wk, double *X, std::complex<double> *A, std::complex<double> *B);

template<class F>
void TanhSinh(const F f0, double a, double b, double eps, std::complex<double> &val);

template<class F>
void TruncIndex(const F f, double eh, int &n, std::complex<double> &s, double sigma, double eta, double a, double b);

template<class F>
void PartSum(const F f, double eh, double e2h, int &n, std::complex<double> &s, double sigma, double eta, double a, double b);

template<class F>
std::complex<double> Term(const F f, double ekh, double sigma, double eta, double a, double b);


// ====== Wrappers for computing the Bessel function ======

void BesselJ_Amos(double v, std::complex<double> x, std::complex<double> a, std::complex<double> &result);

double BesselJ_NextZero(double v, double x);


// ====== Numerical integration wrapper functions via Boost ======

template<class F>
void GaussKronrodBoost(const F f, double a, double b, double tol, std::complex<double> &result);

template<class F>
void TanhSinhBoost(const F f, double a, double b, double tol, std::complex<double> &result);

template<class F>
void ExpSinhBoost(const F f, double a, double tol, std::complex<double> &result);


// ==================================================================================
// Definitions
// ==================================================================================
// Computational drivers
// ==================================================================================

/*! \brief Function to assemble the integrand for a given component of the MGF.*/
inline std::complex<double> SommerfeldIntegrand(SpectralMGF &smgf, double rho, std::complex<double> krho, int component, int order, bool curl)
{

	// Compute spectral MGF or its curl product
	std::complex<double> K;
	smgf.SetRadialWaveNumber(krho);
	if (!curl)
	{		
		smgf.ComputeSpectralMGF();
		K = smgf.GetResult(component);
	}
	else
	{
		smgf.ComputeCurlSpectralMGF();
		K = smgf.GetResultCurl(component);
	}

	// Compute Bessel function
	std::complex<double> Jv;
	BesselJ_Amos(order, krho, rho, Jv);

	return Jv*K*std::pow(krho, order+1)/2.0/M_PI;

}


/*! \brief Function to integrate the spectral MGF for values of krho beyond the spectral near-field, using partition-extrapolation routines.*/
inline std::complex<double> IntegrateSpectralFarField(SpectralMGF &smgf, double rho, int component, int order, double a, bool curl)
{

	double tol = 1.0e-6;
	
	// According to Michalski, Mosig, 2016, the partition-extrapolation procedure should begin at the first zero of the Bessel function that occurs after a - call this point b. A separate finite-domain integration has to be performed to "bridge the gap" between a and b. Pre-tabulated values of Bessel zeros from Abramowitz & Stegun 1972 are used to pick b, available through the Boost library.

	double b = BesselJ_NextZero(order, a*rho)/rho;
	double q = M_PI/rho;

	// MGF integrand functor
	auto f = [&rho, &smgf, &component, &order, &curl](std::complex<double> krho) -> std::complex<double>
		{
			return SommerfeldIntegrand(smgf, rho, krho, component, order, curl);
		};

	// For small rho, Bessel function oscillations are slow, so use quadrature
	// if (rho/(2.0*M_PI/std::real(smgf.lm->k_max)) < 0.01)
	if (rho < 1.0e-15)
	{
		std::complex<double> result = 0.0;
		GaussKronrodBoost(f, a, std::numeric_limits<double>::infinity(), tol, result);
		return result;
	};
	
	std::complex<double> bridge = 0.0;
	if (b > a)
		// TanhSinh(f, a, b, tol, bridge);
		GaussKronrodBoost(f, a, b, tol, bridge);

	std::complex<double> result = 0.0;
	PartExtrap(f, b, q, tol, result);

	
	// DEBUGGING

	// if (curl && component == 2)
	// std::cout << result << std::endl;

	// if (!std::isfinite(std::abs(bridge)))
	// 	std::cout << "bridge: " << bridge << "; rho: " << rho << "; component: " << component << std::endl;
	// if (!std::isfinite(std::abs(result)))
	// 	std::cout << "result: " << result << "; rho: " << rho << "; component: " << component << std::endl;

	// \todo Extremely small distances give overflow / underflow issues - needs to be addressed
	// if (!std::isfinite(std::abs(bridge)))
	// 	bridge = 0.0;
	// if (!std::isfinite(std::abs(result)))
		// result = 0.0;

	return result + bridge;

}


/*! \brief Function to integrate the spectral MGF for values of krho in the spectral near-field, using a rectangular path deformation to avoid singularities.*/
inline std::complex<double> IntegrateSpectralNearField(SpectralMGF &smgf, double rho, int component, int order, double a, bool curl)
{

	// if (std::abs(smgf.lm->k_max) < 1.0e-15)
	// 	return 0.0;

	// Integrate along three paths in the complex plane: first, upwards along the imaginary axis, then right-wards until all poles are left behind, and then back down towards the real axis.

	const std::complex<double> J (0.0, 1.0);
	
	double h = 1.0e-2*std::abs(smgf.lm->k_max);
	// double h = 1.0e-3*std::abs(smgf.lm->k_max);
	// double h = 1.0e-1*std::abs(smgf.lm->k_max);
	double krho_real, krho_imag;
	
	// MGF integrand functor for a fixed imaginary part in the complex plane
	auto f_real = [&rho, &smgf, &component, &order, &curl, &krho_imag, &J](double _krho_real) -> std::complex<double>
		{
			std::complex<double> krho = _krho_real + J*krho_imag;
			return SommerfeldIntegrand(smgf, rho, krho, component, order, curl);
		};

	// MGF integrand functor for a fixed real part in the complex plane
	auto f_imag = [&rho, &smgf, &component, &order, &curl, &krho_real, &J](double _krho_imag) -> std::complex<double>
		{
			std::complex<double> krho = krho_real + J*_krho_imag;
			return SommerfeldIntegrand(smgf, rho, krho, component, order, curl);
		};

	double tol = 1.0e-6;
	double zero_offset = 0.0e-3;
	
	// First path
	std::complex<double> first_path = 0.0;
	krho_real = zero_offset;
	GaussKronrodBoost(f_imag, 0.0, h, tol, first_path);

	// Second path
	std::complex<double> second_path = 0.0;
	krho_imag = h;
	GaussKronrodBoost(f_real, zero_offset, a, tol, second_path);
	// TanhSinhBoost(f_real, zero_offset, a, tol, second_path);
	
	// First path
	std::complex<double> third_path = 0.0;
	krho_real = a;
	GaussKronrodBoost(f_imag, h, 0.0, tol, third_path);
	
	return (first_path + second_path + third_path);

}


// ==================================================================================
// Partition-extrapolation integration routines
// ==================================================================================

/*! \brief Driver function for integrating the spectral far-field (SFF) tail integral of the MGF using techniques from Michalski, Mosig, 2016.*/
template<class F>
void PartExtrap(const F f, double a, double q, double tol, std::complex<double> &val)
{
	
	// Initialize
	int kmax = 10;
	double eps = 1.0e-6;
	std::complex<double> old = 1.0e300, u, w, s = 0.0;
	std::complex<double> A [kmax + 1], B [kmax + 1];
	double X [kmax + 2];
	X[0] = a;

	// Begin extrapolated integration over partitions
	for (int kk = 1; kk <= kmax+1; kk++)
	{
		// Execute quadrature over the current partition
		X[kk] = X[kk-1] + q;
		TanhSinh(f, X[kk-1], X[kk], eps, u);

		if (std::abs(u) == 0.0)
			continue;

		// Execute extrapolation		
		s += u;
		w = u; 
		val = LevinSidi(kk-1, s, w, X, A, B);
		if (kk > 0 && std::abs(val - old) < tol*std::abs(val))
			break;
		
		old = val;
	}
	
	return;

}


/*! \brief Levin transformation function with the Sidi W-algorithm extrapolation.*/
inline std::complex<double> LevinSidi(int k, std::complex<double> Sk, std::complex<double> wk, double *X, std::complex<double> *A, std::complex<double> *B)
{
	
	B[k] = 1.0/wk; 
	A[k] = Sk*B[k];

	for (int j = 1; j <= k; j++)
	{
		double d = 1.0/X[k+1] - 1.0/X[k-j+1]; 
		A[k-j] = (A[k-j+1] - A[k-j])/d;
		B[k-j] = (B[k-j+1] - B[k-j])/d;
	}
	
	return A[0]/B[0];

}


/*! \brief Execute double-exponential tanh-sinh quadrature based on Michalski, Mosig, 2016.*/
template<class F>
void TanhSinh(const F f0, double a, double b, double eps, std::complex<double> &val)
{

	// The algorithm requires the functor to accept two arguments to handle end-point singularities. This should not be an issue for us, so we can side-step this requirement by wrapping the actual functor f0 in a dummy functor.
	auto f = [&f0](double c, double d) -> std::complex<double>
		{
			return f0(c+d);
		};

	double eta = 1.0;
	int maxlev = 5;

	double sigma = (b - a)/2.0;
	double gamma = (b + a)/2.0;
	std::complex<double> s = eta*f(gamma, 0.0);
	
	double h = 0.5;
	double eh = std::exp(h);

	int n;
	TruncIndex(f, eh, n, s, sigma, eta, a, b);
	std::complex<double> old = sigma*h*s;

	for (int m = 1; m <= maxlev; m++)
	{
		double e2h = eh;
		h /= 2.0;
		eh = std::exp(h);
		
		PartSum(f, eh, e2h, n, s, sigma, eta, a, b);
		val = old/2.0 + sigma*h*s;
		
		if (std::abs(val - old) <= eps*std::abs(val))
			break;
		
		old = val;
		n *= 2;
	}

	return;

}


/*! \brief Helper function for tanh-sinh quadrature.*/
template<class F>
void TruncIndex(const F f, double eh, int &n, std::complex<double> &s, double sigma, double eta, double a, double b)
{

	int nmax = 24;
	double kappa = 1.0e-15;
	double ekh = eh;

	for (n = 0; n <= nmax; n++)
	{
		std::complex<double> t = Term(f, ekh, sigma, eta, a, b);
		s += t;

		if (std::abs(t) <= kappa*std::abs(s))
			break;

		ekh *= eh;
	}

	n -= 1;
	
	return;
	
}


/*! \brief Helper function for tanh-sinh quadrature.*/
template<class F>
void PartSum(const F f, double eh, double e2h, int &n, std::complex<double> &s, double sigma, double eta, double a, double b)
{
	
	double ekh = eh;
	s = Term(f, ekh, sigma, eta, a, b);

	for (int k = 2; k <= n; k++)
	{
		ekh *= e2h;
		s += Term(f, ekh, sigma, eta, a, b);
	}

	return;
	
}


/*! \brief Helper function for tanh-sinh quadrature.*/
template<class F>
std::complex<double> Term(const F f, double ekh, double sigma, double eta, double a, double b)
{
	
	double q = std::exp(-eta*(ekh - 1.0/ekh));
	double delta = 2.0*q/(1.0 + q);
	double w = eta*(ekh + 1.0/ekh)*delta/(1.0 + q);
	std::complex<double> t = w*(f(a, sigma*delta) + f(b, -sigma*delta));		

	return t;
	
}


// ==================================================================================
// Wrappers for computing the Bessel function
// ==================================================================================

/*! \brief Compute the Bessel function of first kind, J_v(a*x), using the external library AmosBessel.*/
inline void BesselJ_Amos(double v, std::complex<double> x, std::complex<double> a, std::complex<double> &result)
{

	x *= a;
	__complex__ double x_in;
	__real__ x_in = x.real();
	__imag__ x_in = x.imag();
	
	int nOrders = 1;
	int scale = 0;
	
	int info = AmosBessel('J', x_in, v, nOrders, scale, (__complex__ double *) &result);

	if (info != 0)
		std::cout << "[WARNING] BesselJ_Amos(): AmosBessel() returned " << info << std::endl;
	
	return;

}


/*! \brief Compute the next zero of the Bessel function of first kind, J_v, after a given point x, using the external library Boost.*/
inline double BesselJ_NextZero(double v, double x)
{

	// Use tabulated values, if possible
	if ((int)v == 0)
		return BesselJ0_NextZero(x);
	else if ((int)v == 1)
		return BesselJ1_NextZero(x);

	double result = -1.0;
	int ii = 1;
	
	while (result <= x)
	{
		result = boost::math::cyl_bessel_j_zero((double)v, ii);
		ii++;
	}

	return result;

}


// ==================================================================================
// Numerical integration wrapper functions via Boost
// ==================================================================================

/*! \brief Wrapper function that uses Boost's Gauss-Kronrod quadrature routine to integrate a given complex-valued function along the real line. The user-provided functor should accept a single double and return a complex double.*/
template<class F>
void GaussKronrodBoost(const F f, double a, double b, double tol, std::complex<double> &result)
{
	
	const int order = 31, max_levels = 15;
	double error, multiplier = 1.0;

	// Make sure the integration path is from the smaller number to the larger one, for Boost compatibility.
	if (a > b)
	{
		double c = b;
		b = a;
		a = c;
		multiplier = -1.0;
	}

	typedef typename std::complex<double>::value_type value_type;
	boost::math::quadrature::gauss_kronrod<value_type, order> integrator;
	
	result = integrator.integrate(f, a, b, max_levels, tol, &error);
	result *= multiplier;

	return;
	
}


/*! \brief Wrapper function that uses Boost's double-exponential tanh-sinh quadrature routine to integrate a given complex-valued function along the real line. The user-provided functor should accept a single double and return a complex double.*/
template<class F>
void TanhSinhBoost(const F f, double a, double b, double tol, std::complex<double> &result)
{
	
	double multiplier = 1.0;

	// Make sure the integration path is from the smaller number to the larger one, for Boost compatibility.
	if (a > b)
	{
		double c = b;
		b = a;
		a = c;
		multiplier = -1.0;
	}

	typedef typename std::complex<double>::value_type value_type;
	boost::math::quadrature::tanh_sinh<value_type> integrator;

	result = integrator.integrate(f, a, b, tol);
	result *= multiplier;

	return;
	
}


/*! \brief Wrapper function that uses Boost's double-exponential exp-sinh quadrature routine to integrate a given complex-valued function along the real line, from a given start point to positive infinity. The user-provided functor should accept a single double and return a complex double.*/
template<class F>
void ExpSinhBoost(const F f, double a, double tol, std::complex<double> &result)
{
	
	typedef typename std::complex<double>::value_type value_type;
	boost::math::quadrature::exp_sinh<value_type> integrator;

	result = integrator.integrate(f, a, std::numeric_limits<double>::infinity(), tol);

	return;
	
}


#endif





