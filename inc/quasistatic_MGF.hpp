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

/******************************** quasistatic_MGF.hpp ******************************

 * Routines for computing the spectral and spatial quasistatic approximations to the
 * multilayer Green's function (MGF) based on Michalski & Zheng's Formulation-C.
 *
 * References:
 * Michalski, Zheng, "Electromagnetic Scattering and Radiation by Surfaces of 
 * Arbitrary Shape in Layered Media, Part I: Theory", IEEE-TAP, 1990
 * Simsek, Liu, Wei, "Singularity Subtraction for Evaluation of Green's Functions
 * for Multilayer Media", IEEE-MTT, 2006
 *
 * Author: Shashwat Sharma
 * Created on: Mar 19, 2020

 ***********************************************************************************/


#ifndef QMGF_H
#define QMGF_H


#include <complex>
#include <vector>
#include <array>

#include "layers.hpp"


/*

Usage:

0. Create and initialize spectral MGF generator object
       QuasistaticMGF qmgf;
	   qmgf.Initialize(&layer_manager_object, frequency, <extraction_settings>);

1. Set source and observation layer indices:
       qmgf.SetLayers(src_layer_index, obs_layer_index);

2. Set source and test z points:
       qmgf.SetSourcePoint(z_src);
       qmgf.SetObservationPoint(z_obs);

3. Set lateral wave vector (needed only for spectral quasistatic MGF)
       qmgf.SetRadialWaveNumber(krho);

4. Compute spectral quasistatic MGF components:
       qmgf.ComputeSpectralQMGF();

6. Compute spatial quasistatic MGF components:
       qmgf.ComputeQMGF();

7. Extract output for desired components (spectral or spatial, depending on previous call):
       qmgf.GetResult(component);

*/


/*! \brief Struct for interface Fresnel reflection coefficients.*/
struct Fresnel
{
	std::complex<double> Fe, Fh;
};


/*! \brief Class to generate spectral and spatial quasistatic MGF components.*/
class QuasistaticMGF
{
public:

	// ============ Interface ============

    void Initialize(LayerManager *_lm, double _f, std::vector<bool> _components = {1, 1, 1, 1, 1}, bool _extract_singularities = false);
    void SetLayers(int _i, int _m);
	void SetCurlToGEM();
	void SetCurlToGHJ();

	void ComputeQMGF_Spectral(double z, double zp, std::complex<double> krho);
	std::complex<double> GetResult_Spectral(int component);

	void ComputeCurlQMGF_Spectral(double z, double zp, std::complex<double> krho);
	std::complex<double> GetResultCurl_Spectral(int component);
	
	void ComputeQMGF_Spatial(double z, double zp, double rho);
	std::complex<double> GetResult_Spatial(int component);

	void ComputeCurlQMGF_Spatial(double z, double zp, double rho);
	std::complex<double> GetResultCurl_Spatial(int component);

	void ComputeSingularityFactors(double z, double zp, double rho);
	std::complex<double> GetSingularityFactor(int component);
	std::complex<double> GetSingularityFactorCurl(int component);

	
    // ============ Precomputation ============

	void ComputeFresnelCoefficients();
	Fresnel GetFresnelCoefficient_Upward(int idx);
	Fresnel GetFresnelCoefficient_Downward(int idx);


	// ============ Helpers ============

	std::complex<double> ComputeAxialWaveNumber(std::complex<double> k, std::complex<double> krho);
	void CheckQuadrants(std::complex<double> &kz_test);


	// ============ Computational drivers ============

	void ComputeKii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &K);
	void ComputeKmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &K);

	void ComputeKii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 5> &K);
	void ComputeKmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 5> &K);

	void ComputeFii(double z, double zp, double rho);
	void ComputeFmi(double z, double zp, double rho);
	

	// ============ Computational drivers - curl ============

	void ComputeCurlGEMii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK);
	
	void ComputeCurlGEMmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK);

	void ComputeCurlGEMii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK);
	
	void ComputeCurlGEMmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK);


	// ============ Computational helpers ============

	void ComputeExpTermsKii_kz(double z, double zp, std::complex<double> kzi, std::array<std::complex<double>, 5> &exp_terms);
	void ComputeExpTermsKii_krho(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &exp_terms);
	
	void ComputeExpTermsKmi_krho(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 10> &exp_terms);
	void ComputeExpTermsKmi_kz(double z, double zp, std::complex<double> kzi, std::array<std::complex<double>, 10> &exp_terms);

	void ComputeExpTermsKii_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 5> &exp_terms);
	void ComputeExpTermsKii_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 5> &exp_terms);
	
	void ComputeExpTermsKmi_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 10> &exp_terms);
	void ComputeExpTermsKmi_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 10> &exp_terms);

	void ComputeDistanceTermsKii(double z, double zp, double rho, std::array<double, 5> &ksi);
	void ComputeDistanceTermsKmi(double z, double zp, double rho, std::array<double, 10> &ksi);

	void ComputeGammaTerms(double z, double zp, std::array<double, 5> &gamma);
	void ComputeTauTerms(double z, std::array<double, 2> &tau);
	void ComputeAlphaTerms(double z, double zp, std::array<double, 10> &alpha);

	void ComputeM_Upward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh, double &_nu);
	void ComputeM_Downward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh, double &_nu);

	std::complex<double> ComputeNonsingularHGF(double r, std::complex<double> k);
	bool AxiallyCoincident(double z, double zp, double tol = 1.0e-6);
	bool LaterallyCoincident(double rho, double tol = 1.0e-6);
	

	// ============ Computational helpers - curl ============

	void ComputeTermsCurlKii_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms);
	void ComputeTermsCurlKii_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms);
	void ComputeTermsCurlKii_SpatialS2(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms);
	
	void ComputeTermsCurlKmi_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms);
	void ComputeTermsCurlKmi_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms);
	void ComputeTermsCurlKmi_SpatialS2(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms);

	void SwapSourceAndObserver(double &z, double &zp);
	
	
	// ============ Storage ============

	LayerManager *lm;
    double f, omega;

	std::vector<bool> components = {1, 1, 1, 1, 1};
	std::vector<bool> components_curl = {1, 1, 1, 1};
	bool extract_singularities = false;
	bool curl_GEM = true;
	
	std::vector<Fresnel> Rrs, Rls;

	int i, m;
	std::complex<double> ki, km, epsi, epsm, mui, mum;

	bool initialized = false;
	bool layers_set = false;
	
	std::complex<double> P;
	double nu;
	Fresnel Rm;

	std::complex<double> MVe, MVh;
	std::complex<double> MIe, MIh;

	std::array<std::complex<double>, 5> K_spectral;
	std::array<std::complex<double>, 5> K_spatial;
	std::array<std::complex<double>, 5> F;
	
	std::array<std::complex<double>, 4> curlK_spectral;
	std::array<std::complex<double>, 4> curlK_spatial;
	std::array<std::complex<double>, 4> curlF;

	
};

#endif
