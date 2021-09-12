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

/********************************** spectral_MGF.hpp *******************************

 * Routines for computing the spectral multilayer Green's function (MGF) based on
 * Michalski & Zheng's Formulation-C.
 *
 * References:
 * Michalski, Zheng, "Electromagnetic Scattering and Radiation by Surfaces of 
 * Arbitrary Shape in Layered Media, Part I: Theory", IEEE-TAP, 1990
 *
 * Author: Shashwat Sharma
 * Created on: Feb 15, 2019

 ***********************************************************************************/


#ifndef SPECTRALMGF_H
#define SPECTRALMGF_H


#include <complex>
#include <vector>
#include <array>

#include "layers.hpp"
#include "quasistatic_MGF.hpp"


/*

Usage:

0. Create and initialize spectral MGF generator object
       SpectralMGF smgf;
	   smgf.Initialize(&layer_manager_object, frequency, <extraction_settings>);

1. Set source and observation layer indices:
       smgf.SetLayers(src_layer_index, obs_layer_index);

2. Set source and test z points:
       smgf.SetSourcePoint(z_src);
       smgf.SetObservationPoint(z_obs);

3. Set lateral wave vector
       smgf.SetRadialWaveNumber(krho);

4. Compute spectral MGF components:
       smgf.ComputeSpectralMGF();

5. Extract output for desired components:
       smgf.GetResult(component);

Notes:
  
- If possible, loop first over krho, and then over z points, so that the same krho precomputations (which are significant) can be reused as much as possible.

- If the above is possible, then make sure to call SetSourcePoint() and SetObservationPoint() in the outermost possible loops, because that would allow additional precomputations to be done. For example, if first looping over source points and then over observation points, then call SetSourcePoint() outside the inner observation point loop, not inside it. That way, all source-related terms will only be computed once for all observations points. Note that this only applies when source and observation points are in different layers, and if krho is set before setting the source and observation points.

- Importantly, the layer indices must be set before the radial wave number and the source and observation points, because all precomputations rely on knowledge of source and observation layers.

*/


/*! \brief Struct for layer and interface impedance.*/
struct Impedance
{
	std::complex<double> Ze, Ye, Zh, Yh;
};


/*! \brief Struct for directional layer interface reflection coefficients.*/
struct Gamma
{
	std::complex<double> Ge, Gh;
};


/*! \brief Class to generate spectral MGF components.*/
class SpectralMGF
{
public:

	// ============ Interface ============

    void Initialize(LayerManager *_lm, double _f, std::vector<bool> _components = {1, 1, 1, 1, 1}, bool _extract_quasistatic = false, bool _extract_homogeneous = false);
    void SetLayers(int _i, int _m);
	void SetSourcePoint(double _zp);
	void SetObservationPoint(double _z);
	void SetRadialWaveNumber(std::complex<double> _krho);
	void SetCurlToGEM();
	void SetCurlToGHJ();
	
	void ComputeSpectralMGF();
	std::complex<double> GetResult(int component);

	void ComputeCurlSpectralMGF();
	std::complex<double> GetResultCurl(int component);

	
    // ============ Precomputation ============
	
    void ComputeHalfspaceImpedance(std::complex<double> krho);
	void ComputeLayerImpedance(std::complex<double> krho, int idx, Impedance &Z, std::complex<double> *_kz = NULL);
	
    void ComputeInterfaceImpedance_Upward(std::complex<double> krho, int idx, Impedance &Zr);
	void ComputeInterfaceImpedance_Downward(std::complex<double> krho, int idx, Impedance &Zl);

	void ComputeReflectionCoefficient_Upward(std::complex<double> krho, int idx, Gamma &Gr);
	void ComputeReflectionCoefficient_Downward(std::complex<double> krho, int idx, Gamma &Gl);


	// ============ Helpers ============

	std::complex<double> ComputeAxialWaveNumber(std::complex<double> k, std::complex<double> krho);
	void CheckQuadrants(std::complex<double> &kz_test);


	// ============ Computational drivers ============

	void ComputeKii(std::array<std::complex<double>, 5> &K);
	void ComputeKmi(std::array<std::complex<double>, 5> &K);


	// ============ Computational helpers ============

	void ComputeGVii(double _z, double _zp, std::complex<double> &GVe, std::complex<double> &GVh);
	void ComputeGIii(double _z, double _zp, std::complex<double> &GIe, std::complex<double> &GIh);
	void ComputeWVii(double _z, double _zp, std::complex<double> &WVe, std::complex<double> &WVh);
	void ComputeWIii(double _z, double _zp, std::complex<double> &WIe, std::complex<double> &WIh);

	void ComputeExpTermsGii(double _z, double _zp, std::array<std::complex<double>, 3> &_exp_terms, std::complex<double> &_cos_term, std::complex<double> &_sin_term);

	void ComputeTV(double _z, std::complex<double> &_TVe, std::complex<double> &_TVh);
	void ComputeTI(double _z, std::complex<double> &_TIe, std::complex<double> &_TIh);
	void ComputeTrfTerms(double _z, std::array<std::complex<double>, 2> &_trf_terms);

	void ComputeM_Upward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh);
	void ComputeM_Downward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh);


	// ============ Computational drivers - curl ============

	void ComputeCurlGEMii(std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJii(std::array<std::complex<double>, 4> &curlK);
	
	void ComputeCurlGEMmi(std::array<std::complex<double>, 4> &curlK);
	void ComputeCurlGHJmi(std::array<std::complex<double>, 4> &curlK);


	// ============ Computational helpers - curl ============
	
	void ComputeVVhii(double _z, double _zp, std::complex<double> &_VVh);
	void ComputeVVeii(double _z, double _zp, std::complex<double> &_VVe);

	void ComputeIIhii(double _z, double _zp, std::complex<double> &_IIh);
	void ComputeIIeii(double _z, double _zp, std::complex<double> &_IIe);

	void ComputeVIhii(double _z, double _zp, std::complex<double> &_VIh);
	void ComputeVIeii(double _z, double _zp, std::complex<double> &_VIe);
	void ComputeIVeii(double _z, double _zp, std::complex<double> &_IVe);

	void ComputeExpTermsGii_dz(double _z, double _zp, std::array<std::complex<double>, 3> &_exp_terms, std::complex<double> &_cos_term, std::complex<double> &_sin_term);

	void SwapSourceAndObserver();
	
	
	// ============ Storage ============

	LayerManager *lm;
    double f, omega;

	QuasistaticMGF qmgf;

	std::vector<bool> components = {1, 1, 1, 1, 1};
	std::vector<bool> components_curl = {1, 1, 1, 1};
	bool extract_quasistatic = false;
	bool extract_homogeneous = false;
	
    int i, m;
    double zp, z;
	std::complex<double> ki, km, epsi, epsm, mui, mum;
	
	std::array<std::complex<double>, 3> exp_terms;
	std::array<std::complex<double>, 2> trf_terms;
	std::complex<double> cos_term, sin_term;

	bool initialized = false;
	bool layers_set = false;
	bool src_point_set = false;
	bool obs_point_set = false;
	bool precomputations_done = false;
	bool exp_terms_computed = false;
	bool src_terms_computed = false;
	bool obs_terms_computed = false;
	bool trf_terms_computed = false;

	std::complex<double> krho, krhosq;
	std::complex<double> kz0, kzi, kzm, jkzi, jkzm;

	Impedance Z_top, Z_bot;
	Impedance Zi, Zm, Zdi;
	Gamma Gri, Gli, Gm;
	
	std::complex<double> De, Dh, P;

	std::complex<double> MVe, MVh;
	std::complex<double> MIe, MIh;

	std::complex<double> GVe, GVh;
	std::complex<double> GIe, GIh;
	
	std::complex<double> WVe, WVh;
	std::complex<double> WIe, WIh;
	
	std::complex<double> TVe, TVh;
	std::complex<double> TIe, TIh;

	std::array<std::complex<double>, 5> K;

	// ------ Curl ------

	bool curl_GEM = true;
	
	std::complex<double> VVe, VVh;
	std::complex<double> IIe, IIh;
	std::complex<double> IVe, VIh;
	
	std::array<std::complex<double>, 4> curlK;
	

};

#endif
