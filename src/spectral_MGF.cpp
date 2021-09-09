/********************************** spectral_MGF.cpp *******************************

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


#include <cmath>
#include <complex>
#include <stdexcept>

#include "spectral_MGF.hpp"
#include "constants.hpp"

using namespace strata;


// ==================================================================================
// Interface
// ==================================================================================

/*! \brief Initialize the spectral MGF engine for a given frequency and layered environment.*/
void SpectralMGF::Initialize(LayerManager *_lm, double _f, std::vector<bool> _components, bool _extract_quasistatic, bool _extract_homogeneous)
{
	
	lm = _lm;
	f = _f;
	omega = 2.0*M_PI*f;

	components = _components;
	extract_quasistatic = _extract_quasistatic;
	extract_homogeneous = _extract_homogeneous;
	
	if (!lm->layers_processed)
		lm->ProcessLayers(f);

	if (extract_quasistatic)
		qmgf.Initialize(lm, f, components);

	initialized = true;
	layers_set = false;
	src_point_set = false;
	obs_point_set = false;
	precomputations_done = false;
	exp_terms_computed = false;
	trf_terms_computed = false;
	
	return;
	
}


/*! \brief Set source and observation layers and execute related precomputations.*/
void SpectralMGF::SetLayers(int _i, int _m)
{

	// ------ Sanity checks ------
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] SpectralMGF::SetLayers(): Must call Initialize() before setting layer indices.");
	}

	// if (_i < -1 || _i > lm->layers.size())
	// {
	// 	throw std::domain_error("[ERROR] QuasistaticMGF::SetLayers(): Source layer index is out of bounds.");
	// }

	// if (_m < -1 || _m > lm->layers.size())
	// {
	// 	throw std::domain_error("[ERROR] QuasistaticMGF::SetLayers(): Observation layer index is out of bounds.");
	// }

	// ---------------------------
	
	i = _i;
	m = _m;

	if (i == -1)
	{
		ki = lm->k_top;
		epsi = lm->eps_top;
		mui = lm->mu_top;
	}
	else if (i == lm->layers.size())
	{
		ki = lm->k_bot;
		epsi = lm->eps_bot;
		mui = lm->mu_bot;
	}
	else
	{
		ki = lm->k[i];
		epsi = lm->eps[i];
		mui = lm->mu[i];
	}

	if (m == -1)
	{
		km = lm->k_top;
		epsm = lm->eps_top;
		mum = lm->mu_top;
	}
	else if (m == lm->layers.size())
	{
		km = lm->k_bot;
		epsm = lm->eps_bot;
		mum = lm->mu_bot;
	}
	else
	{
		km = lm->k[m];
		epsm = lm->eps[m];
		mum = lm->mu[m];
	}

	if (extract_quasistatic)
		qmgf.SetLayers(i, m);

	layers_set = true;
	src_point_set = false;
	obs_point_set = false;
	precomputations_done = false;
	exp_terms_computed = false;
	trf_terms_computed = false;
	
	return;
	
}


/*! \brief Set the curl of the MGF to G_EM, for EFIE-type formulations.*/
void SpectralMGF::SetCurlToGEM()
{
	curl_GEM = true;
	if (extract_quasistatic)
		qmgf.SetCurlToGEM();
	return;
}


/*! \brief Set the curl of the MGF to G_HJ, for MFIE-type formulations.*/
void SpectralMGF::SetCurlToGHJ()
{
	curl_GEM = false;
	if (extract_quasistatic)
		qmgf.SetCurlToGHJ();
	return;
}


/*! \brief Set source z-value.*/
void SpectralMGF::SetSourcePoint(double _zp)
{

	zp = _zp;

	// Precompute source-dependent terms, if possible
	src_terms_computed = false;
	
	if (precomputations_done && i != m)
	{
		double _z;
		
		if (m < i)
			// _z = lm->layers[i].zmax;
			_z = lm->GetZmax(i);
		else if (m > i)
			// _z = lm->layers[i].zmin;
			_z = lm->GetZmin(i);

		if (components[0] || components[2] || components[4])
			ComputeGVii(_z, zp, GVe, GVh);
		if (components[1] || components[3])
			ComputeGIii(_z, zp, GIe, GIh);
		
		src_terms_computed = true;
	}
	
	src_point_set = true;
	exp_terms_computed = false;

	return;
	
}


/*! \brief Set observation z-value.*/
void SpectralMGF::SetObservationPoint(double _z)
{
	
	z = _z;

	// Precompute observer-dependent terms, if possible
	obs_terms_computed = false;
	
	if (precomputations_done && i != m)
	{
		if (components[0] || components[1] || components[4])
			ComputeTV(z, TVe, TVh);
		if (components[2] || components[3])
			ComputeTI(z, TIe, TIh);
		
		obs_terms_computed = true;
	}
	
	obs_point_set = true;
	exp_terms_computed = false;
	trf_terms_computed = false;

	return;
	
}


/*! \brief Interface function to set the radial (lateral) wave number, krho, and execute relevant precomputations.*/
void SpectralMGF::SetRadialWaveNumber(std::complex<double> _krho)
{

	// ------ Check that initializations have been done ------
	
	if (!layers_set)
	{
		throw std::logic_error("[ERROR] SpectralMGF::SetRadialWaveNumber(): Must call SetLayers() before setting the radial wave number.");
	}

	// -------------------------------------------------------

	krho = _krho;
	krhosq = krho*krho;

	kz0 = ComputeAxialWaveNumber(lm->k_min, krho);
	kzi = ComputeAxialWaveNumber(ki, krho);
	kzm = ComputeAxialWaveNumber(km, krho);
	jkzi = J*kzi;
	jkzm = J*kzm;

	ComputeHalfspaceImpedance(krho);
	
	ComputeLayerImpedance(krho, i, Zi, &kzi);
	ComputeLayerImpedance(krho, m, Zm, &kzm);

	ComputeReflectionCoefficient_Upward(krho, i-1, Gri);
	ComputeReflectionCoefficient_Downward(krho, i, Gli);

	// std::complex<double> exp_term = std::exp(-jkzi*2.0*lm->layers[i].h);
	std::complex<double> exp_term = std::exp(-jkzi*2.0*lm->GetHeight(i));

	if (i == -1 || i == (int)lm->layers.size())
		exp_term = 0.0;

	De = 1.0 - Gli.Ge*Gri.Ge*exp_term;
	Dh = 1.0 - Gli.Gh*Gri.Gh*exp_term;

	if (m < i)
	{
		// Observation layer is above the source layer
		ComputeInterfaceImpedance_Upward(krho, i-1, Zdi);
		ComputeReflectionCoefficient_Upward(krho, m-1, Gm);
		ComputeM_Upward(MVe, MVh, MIe, MIh);
	}
	else if (m > i)
	{		
		// Observation layer is below the source layer
		ComputeInterfaceImpedance_Downward(krho, i, Zdi);
		ComputeReflectionCoefficient_Downward(krho, m, Gm);
		ComputeM_Downward(MVe, MVh, MIe, MIh);
	}

	precomputations_done = true;
	exp_terms_computed = false;
	trf_terms_computed = false;
	src_terms_computed = false;
	obs_terms_computed = false;
	
	return;
	
}


/*! \brief Driver function to compute desired components of the spectral MGF.*/
void SpectralMGF::ComputeSpectralMGF()
{

	// ------ Check that precomputations have been done ------
	
	if (!src_point_set || !obs_point_set)
	{
		throw std::logic_error("[ERROR] SpectralMGF::ComputeSpectralMGF(): Must call SetSourcePoint() and SetObservationPoint() before computing the spectral MGF.");
	}
	
	if (!precomputations_done)
	{
		throw std::logic_error("[ERROR] SpectralMGF::ComputeSpectralMGF(): Must call SetRadialWaveNumber() before computing the spectral MGF.");
	}
	
	// -------------------------------------------------------

	std::fill(K.begin(), K.end(), 0.0);
	
	if (i == m)
		ComputeKii(K);
	else
		ComputeKmi(K);

	K[0] *= 1.0/J/omega/mu0;
	K[1] *= -mui/J/omega/epsm/mu0;
	K[2] *= -1.0/J/omega/mu0;
	K[3] *= mum/J/omega/epsi/mu0;
	K[4] *= J*omega*eps0;

	if (extract_quasistatic)
	{
		qmgf.ComputeQMGF_Spectral(z, zp, krho);

		K[0] -= qmgf.GetResult_Spectral(0);
		K[1] -= qmgf.GetResult_Spectral(1);
		K[2] -= qmgf.GetResult_Spectral(2);
		K[3] -= qmgf.GetResult_Spectral(3);
		K[4] -= qmgf.GetResult_Spectral(4);
	}
	else if (extract_homogeneous)
	{
		K[0] -= (mui/mu0)*std::exp(-jkzi*std::abs(z - zp))/(2.0*jkzi);
		K[3] -= (mui/mu0)*std::exp(-jkzi*std::abs(z - zp))/(2.0*jkzi);
		K[4] -= (eps0/epsi)*std::exp(-jkzi*std::abs(z - zp))/(2.0*jkzi);
	}
	
	return;
	
}


/*! \brief Extract computed results.*/
std::complex<double> SpectralMGF::GetResult(int component)
{
	return K[component];
}


/*! \brief Driver function to compute desired components of the curl of the spectral MGF.*/
void SpectralMGF::ComputeCurlSpectralMGF()
{

	// ------ Check that precomputations have been done ------
	
	if (!src_point_set || !obs_point_set)
	{
		throw std::logic_error("[ERROR] SpectralMGF::ComputeCurlSpectralMGF(): Must call SetSourcePoint() and SetObservationPoint() before computing the spectral MGF.");
	}
	
	if (!precomputations_done)
	{
		throw std::logic_error("[ERROR] SpectralMGF::ComputeCurlSpectralMGF(): Must call SetRadialWaveNumber() before computing the spectral MGF.");
	}
	
	// -------------------------------------------------------

	if (z > zp)
		P = 1.0;
	else if (z < zp)
		P = -1.0;
	else
		P = 0.0;

	std::fill(curlK.begin(), curlK.end(), 0.0);

	if (curl_GEM)
	{
		if (i == m)
			ComputeCurlGEMii(curlK);
		else
			ComputeCurlGEMmi(curlK);
	}
	else
	{
		if (i == m)
			ComputeCurlGHJii(curlK);
		else
			ComputeCurlGHJmi(curlK);
	}

	if (extract_quasistatic)
	{
		qmgf.ComputeCurlQMGF_Spectral(z, zp, krho);

		curlK[0] -= qmgf.GetResultCurl_Spectral(0);
		curlK[1] -= qmgf.GetResultCurl_Spectral(1);
		curlK[2] -= qmgf.GetResultCurl_Spectral(2);
		curlK[3] -= qmgf.GetResultCurl_Spectral(3);
	}
	
	return;
	
}


/*! \brief Extract computed results.*/
std::complex<double> SpectralMGF::GetResultCurl(int component)
{
	return curlK[component];
}


// ==================================================================================
// Precomputation
// ==================================================================================

/*! \brief Compute top and bottom half-space impedance.*/
void SpectralMGF::ComputeHalfspaceImpedance(std::complex<double> krho)
{

	std::complex<double> kz_top = ComputeAxialWaveNumber(lm->k_top, krho);

	if (!lm->isPEC_top)
	{
		Z_top.Ze = kz_top/(omega*lm->eps_top);
		Z_top.Ye = 1.0/Z_top.Ze;
		Z_top.Zh = (omega*lm->mu_top)/kz_top;
		Z_top.Yh = 1.0/Z_top.Zh;
	}
	else
	{
		Z_top.Ze = 0.0;
		Z_top.Ye = 1.0/Z_top.Ze;
		Z_top.Zh = 0.0;
		Z_top.Yh = 1.0/Z_top.Zh;
	}

	std::complex<double> kz_bot = ComputeAxialWaveNumber(lm->k_bot, krho);

	if (!lm->isPEC_bot)
	{
		Z_bot.Ze = kz_bot/(omega*lm->eps_bot);
		Z_bot.Ye = 1.0/Z_bot.Ze;
		Z_bot.Zh = (omega*lm->mu_bot)/kz_bot;
		Z_bot.Yh = 1.0/Z_bot.Zh;
	}
	else
	{
		Z_bot.Ze = 0.0;
		Z_bot.Ye = 1.0/Z_bot.Ze;
		Z_bot.Zh = 0.0;
		Z_bot.Yh = 1.0/Z_bot.Zh;
	}

	return;
	
}


/*! \brief Compute the partial impedance (Ze, Ye, Zh, Yh) of a given layer. Computes kz, if not already provided (or provided as NULL).*/
void SpectralMGF::ComputeLayerImpedance(std::complex<double> krho, int idx, Impedance &Z, std::complex<double> *_kz)
{

	if (idx < -1 || idx > (int)(lm->layers.size()))
		throw std::out_of_range("[ERROR] SpectralMGF::ComputeLayerImpedance(): Layer index out of bounds; received idx " + std::to_string(idx));

	if (idx == -1)
	{
		ComputeHalfspaceImpedance(krho);
		Z = Z_top;
		return;
	}
	else if (idx == (int)(lm->layers.size()))
	{
		ComputeHalfspaceImpedance(krho);
		Z = Z_bot;
		return;
	}

	std::complex<double> k = lm->k[idx];
	std::complex<double> kz;

	if (_kz == NULL)
		kz = ComputeAxialWaveNumber(k, krho);
	else
		kz = *_kz;
			
	Z.Ze = kz/(omega*lm->eps[idx]);
	Z.Ye = 1.0/Z.Ze;
	Z.Zh = (omega*lm->mu[idx])/kz;
	Z.Yh = 1.0/Z.Zh;

	return;

}


/*! \brief Compute the total impedance (Ze, Ye, Zh, Yh) for a given layer interface, looking upwards. ComputeHalfspaceImpedances() must have been called at some point before this function. By convention, indexing of the interfaces is such that the i^{th} interface is the lower boundary of the i^{th} layer, i.e. the i^{th} layer's bottom surface is the i^{th} interface, and its top surface is the (i-1)^{th} interface. The top surface of the top-most layer is indexed with -1. This is a bit weird, but is done for consistency with Michalski, Zheng, TAP, 1990.*/
void SpectralMGF::ComputeInterfaceImpedance_Upward(std::complex<double> krho, int idx, Impedance &Zr)
{

	if (idx < -1 || idx > (int)(lm->layers.size() - 1))
	{
		throw std::out_of_range("[ERROR] SpectralMGF::ComputeInterfaceImpedance_Upward(): Interface index out of bounds; received idx " + std::to_string(idx));
	}
	
	// Compute kz and thickness factor
	std::complex<double> kz = 0.0, t = 0.0;
	if (idx > -1)
	{
		kz = ComputeAxialWaveNumber(lm->k[idx], krho);
		// t = std::tan(kz*lm->layers[idx].h);
		t = std::tan(kz*lm->GetHeight(idx));
		// t = J*(std::exp(-2.0*J*kz*(lm->layers[idx].h)) - 1.0)/(std::exp(-2.0*J*kz*(lm->layers[idx].h)) + 1.0);
	}

	if (idx == -1)
	{
		// Top-most layer interface
		Zr = Z_top;
	}
	else if (idx == 0 && lm->isPEC_top)
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx, Z, &kz);
		
		Zr.Ze = Z.Ze*J*t;
		Zr.Zh = Z.Zh*J*t;
		Zr.Ye = Z.Ye/(J*t);
		Zr.Yh = Z.Yh/(J*t);
	}
	else
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx, Z, &kz);
		
		// Recursively compute impedances of interfaces above this one
		Impedance Zr_nextabove;
		ComputeInterfaceImpedance_Upward(krho, idx-1, Zr_nextabove);

		// If the tangent is infinity, use L'Hospital's rule, otherwise proceed normally
		if (std::isfinite(std::abs(t)))
		{
			Zr.Ze = Z.Ze*(Zr_nextabove.Ze + J*Z.Ze*t) / (Z.Ze + J*Zr_nextabove.Ze*t);
			Zr.Zh = Z.Zh*(Zr_nextabove.Zh + J*Z.Zh*t) / (Z.Zh + J*Zr_nextabove.Zh*t);
			Zr.Ye = Z.Ye*(Zr_nextabove.Ye + J*Z.Ye*t) / (Z.Ye + J*Zr_nextabove.Ye*t);
			Zr.Yh = Z.Yh*(Zr_nextabove.Yh + J*Z.Yh*t) / (Z.Yh + J*Zr_nextabove.Yh*t);
		}
		else
		{
			Zr.Ze = std::pow(Z.Ze, 2) / Zr_nextabove.Ze;
			Zr.Zh = std::pow(Z.Zh, 2) / Zr_nextabove.Zh;
			Zr.Ye = std::pow(Z.Ye, 2) / Zr_nextabove.Ye;
			Zr.Yh = std::pow(Z.Yh, 2) / Zr_nextabove.Yh;
		}
	}

	return;
	
}


/*! \brief Compute the total impedance (Ze, Ye, Zh, Yh) for a given layer interface, looking downwards. ComputeHalfspaceImpedances() must have been called at some point before this function. By convention, indexing of the interfaces is such that the i^{th} interface is the lower boundary of the i^{th} layer, i.e. the i^{th} layer's bottom surface is the i^{th} interface, and its top surface is the (i-1)^{th} interface. The top surface of the top-most layer is indexed with -1. This is a bit weird, but is done for consistency with Michalski, Zheng, TAP, 1990.*/
void SpectralMGF::ComputeInterfaceImpedance_Downward(std::complex<double> krho, int idx, Impedance &Zl)
{

	if (idx < -1 || idx > (int)(lm->layers.size() - 1))
	{
		throw std::out_of_range("[ERROR] SpectralMGF::ComputeInterfaceImpedance_Downward(): Interface index out of bounds; received idx " + std::to_string(idx));
	}

	// Compute kz and thickness factor
	std::complex<double> kz = 0.0, t = 0.0;
	if (idx < (int)(lm->layers.size() - 1))
	{
		kz = ComputeAxialWaveNumber(lm->k[idx+1], krho);
		// t = std::tan(kz*lm->layers[idx+1].h);
		t = std::tan(kz*lm->GetHeight(idx+1));
		// t = J*(std::exp(-2.0*J*kz*(lm->layers[idx+1].h)) - 1.0)/(std::exp(-2.0*J*kz*(lm->layers[idx+1].h)) + 1.0);
	}
	
	if (idx == (int)(lm->layers.size() - 1))
	{
		// Bottom-most layer interface
	    Zl = Z_bot;
	}
	else if (idx == (int)(lm->layers.size() - 2) && lm->isPEC_bot)
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx+1, Z, &kz);
		
		Zl.Ze = Z.Ze*J*t;
		Zl.Zh = Z.Zh*J*t;
		Zl.Ye = Z.Ye/(J*t);
		Zl.Yh = Z.Yh/(J*t);
	}
	else
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx+1, Z, &kz);
		
		// Recursively compute impedances of interfaces below this one
		Impedance Zl_nextbelow;
		ComputeInterfaceImpedance_Downward(krho, idx+1, Zl_nextbelow);

		// If the tangent is infinity, use L'Hospital's rule, otherwise proceed normally
		if (std::isfinite(std::abs(t)))
		{
			Zl.Ze = Z.Ze*(Zl_nextbelow.Ze + J*Z.Ze*t) / (Z.Ze + J*Zl_nextbelow.Ze*t);
			Zl.Zh = Z.Zh*(Zl_nextbelow.Zh + J*Z.Zh*t) / (Z.Zh + J*Zl_nextbelow.Zh*t);
			Zl.Ye = Z.Ye*(Zl_nextbelow.Ye + J*Z.Ye*t) / (Z.Ye + J*Zl_nextbelow.Ye*t);
			Zl.Yh = Z.Yh*(Zl_nextbelow.Yh + J*Z.Yh*t) / (Z.Yh + J*Zl_nextbelow.Yh*t);
		}
		else
		{
			Zl.Ze = std::pow(Z.Ze, 2) / Zl_nextbelow.Ze;
			Zl.Zh = std::pow(Z.Zh, 2) / Zl_nextbelow.Zh;
			Zl.Ye = std::pow(Z.Ye, 2) / Zl_nextbelow.Ye;
			Zl.Yh = std::pow(Z.Yh, 2) / Zl_nextbelow.Yh;
		}
	}

	return;
	
}


/*! \brief Compute the upward reflection coefficient for a given layer interface.*/
void SpectralMGF::ComputeReflectionCoefficient_Upward(std::complex<double> krho, int idx, Gamma &Gr)
{

	if (idx < -2 || idx > (int)(lm->layers.size() - 1))
	{
		throw std::out_of_range("[ERROR] SpectralMGF::ComputeReflectionCoefficient_Upward(): Layer index out of bounds; received idx " + std::to_string(idx));
	}
	
	if (idx == -2)
	{
		Gr.Ge = 0.0;
		Gr.Gh = 0.0;
	}
	else if (idx == -1 && lm->isPEC_top)
	{
		Gr.Ge = -1.0;
		Gr.Gh = -1.0;
	}
	else
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx+1, Z);
		
		Impedance Zr;
		ComputeInterfaceImpedance_Upward(krho, idx, Zr);
		
		Gr.Ge = (Zr.Ze - Z.Ze)/(Zr.Ze + Z.Ze);
		Gr.Gh = (Zr.Zh - Z.Zh)/(Zr.Zh + Z.Zh);
	}

	return;

}


/*! \brief Compute the downward reflection coefficient for a given layer interface.*/
void SpectralMGF::ComputeReflectionCoefficient_Downward(std::complex<double> krho, int idx, Gamma &Gl)
{

	if (idx < -1 || idx > (int)(lm->layers.size()))
	{
		throw std::out_of_range("[ERROR] SpectralMGF::ComputeReflectionCoefficient_Downward(): Layer index out of bounds; received idx " + std::to_string(idx));
	}
	
	if (idx == (int)(lm->layers.size()))
	{
		Gl.Ge = 0.0;
		Gl.Gh = 0.0;
	}
	else if (idx == (int)(lm->layers.size() - 1) && lm->isPEC_bot)
	{
		Gl.Ge = -1.0;
		Gl.Gh = -1.0;
	}
	else
	{
		Impedance Z;
		ComputeLayerImpedance(krho, idx, Z);

		Impedance Zl;
		ComputeInterfaceImpedance_Downward(krho, idx, Zl);
		
		Gl.Ge = (Zl.Ze - Z.Ze)/(Zl.Ze + Z.Ze);
		Gl.Gh = (Zl.Zh - Z.Zh)/(Zl.Zh + Z.Zh);
	}
	
	return;
	
}


// ==================================================================================
// Helpers
// ==================================================================================


/*! \brief Function to compute kz from k and krho, and make sure it lies in the lower-right complex quadrant.*/
std::complex<double> SpectralMGF::ComputeAxialWaveNumber(std::complex<double> k, std::complex<double> krho)
{
	std::complex<double> kz = std::sqrt(k*k - krho*krho);
	CheckQuadrants(kz);
	return kz;
}


/*! \brief Function to ensure that kz lies in the bottom right quadrant of the complex plane.*/
void SpectralMGF::CheckQuadrants(std::complex<double> &kz_test)
{
	
	if (kz_test.imag() > 0)
		kz_test.imag(-std::abs(std::imag(kz_test)));

	if (kz_test.real() < 0)
		kz_test.real(std::abs(std::real(kz_test)));

	return;

}


// ==================================================================================
// Computational drivers
// ==================================================================================

/*! \brief Compute the spectral MGF when the source and observation layers are the same.*/
void SpectralMGF::ComputeKii(std::array<std::complex<double>, 5> &K)
{

	// ====== Precompute terms ======
	
	if (components[0] || components[4])
		ComputeGVii(z, zp, GVe, GVh);
	
	
	// ====== Assemble required component(s) ======

	if (components[0])
	{
		K[0] = GVh;
	}

	if (components[1])
	{
		ComputeWVii(z, zp, WVe, WVh);
		K[1] = (WVe - (km*km/(kzm*kzm))*WVh)/krhosq;
	}

	if (components[2])
	{
		ComputeWIii(z, zp, WIe, WIh);
		K[2] = ((km*km/(kzm*kzm))*WIe - WIh)/krhosq;
	}

	if (components[3])
	{
		ComputeGIii(z, zp, GIe, GIh);
		K[3] = (GIe - (ki*ki/krhosq)*((kzm*kzm/(km*km))*GIe - GIh));
	}

	if (components[4])
	{
		K[4] = (GVe - GVh)/krhosq;
	}
		
	return;

}


/*! \brief Compute the spectral MGF when the source and observation layers are not the same.*/
void SpectralMGF::ComputeKmi(std::array<std::complex<double>, 5> &K)
{

	// ====== Precompute terms ======

	double _z, _I;
	
	if (m < i)
	{
		// _z = lm->layers[i].zmax;
		_z = lm->GetZmax(i);
		_I = 1.0;
	}
	else if (m > i)
	{
		// _z = lm->layers[i].zmin;
		_z = lm->GetZmin(i);
		_I = -1.0;
	}

	std::complex<double> GVe_ii, GVh_ii, GIe_ii, GIh_ii;

	if (!src_terms_computed)
	{
		if (components[0] || components[2] || components[4])
			ComputeGVii(_z, zp, GVe_ii, GVh_ii);
		if (components[1] || components[3])
			ComputeGIii(_z, zp, GIe_ii, GIh_ii);
		
		src_terms_computed = true;
	}

	if (!obs_terms_computed)
	{
		if (components[0] || components[1] || components[4])
			ComputeTV(z, TVe, TVh);
		if (components[2] || components[3])
			ComputeTI(z, TIe, TIh);
		
		obs_terms_computed = true;
	}


	// ====== Assemble required component(s) ======

	if (components[0])
	{
		std::complex<double> GVh_mi = GVh_ii*TVh;
		K[0] = GVh_mi;
	}

	if (components[1])
	{
		std::complex<double> IVe_mi = _I*Zdi.Ze*GIe_ii*TVe;
		std::complex<double> IVh_mi = _I*Zdi.Zh*GIh_ii*TVh;
		std::complex<double> WVe_mi = (-jkzm/Zm.Ze)*IVe_mi;
		std::complex<double> WVh_mi = (-jkzm/Zm.Zh)*IVh_mi;
		
		K[1] = (WVe_mi - (km*km/(kzm*kzm))*WVh_mi)/krhosq;
	}

	if (components[2])
	{
		std::complex<double> IIe_mi = _I*Zdi.Ye*GVe_ii*TIe;
		std::complex<double> IIh_mi = _I*Zdi.Yh*GVh_ii*TIh;
		std::complex<double> WIe_mi = (-jkzm/Zm.Ye)*IIe_mi;
		std::complex<double> WIh_mi = (-jkzm/Zm.Yh)*IIh_mi;
		
		K[2] = ((km*km/(kzm*kzm))*WIe_mi - WIh_mi)/krhosq;
	}

	if (components[3])
	{
		std::complex<double> GIe_mi = GIe_ii*TIe;
		std::complex<double> GIh_mi = GIh_ii*TIh;
		
		K[3] = (GIe_mi - (ki*ki/krhosq)*((kzm*kzm/(km*km))*GIe_mi - GIh_mi));
	}

	if (components[4])
	{
		std::complex<double> GVe_mi = GVe_ii*TVe;
		std::complex<double> GVh_mi = GVh_ii*TVh;

		K[4] = (GVe_mi - GVh_mi)/krhosq;
	}
	
	return;

}


// ==================================================================================
// Computational helpers
// ==================================================================================

/*! \brief Compute the voltage at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeGVii(double _z, double _zp, std::complex<double> &_GVe, std::complex<double> &_GVh)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_GVe = (Zi.Ze/2.0)*(exp_terms[0] + (1.0/De)*( Gli.Ge*exp_terms[1] + Gri.Ge*exp_terms[2] + 2.0*Gli.Ge*Gri.Ge*cos_term ));

	_GVh = (Zi.Zh/2.0)*(exp_terms[0] + (1.0/Dh)*( Gli.Gh*exp_terms[1] + Gri.Gh*exp_terms[2] + 2.0*Gli.Gh*Gri.Gh*cos_term ));

	return;

}


/*! \brief Compute the current at the observation point z (in layer m) due to a unit-strength voltage source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeGIii(double _z, double _zp, std::complex<double> &_GIe, std::complex<double> &_GIh)
{

	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_GIe = (Zi.Ye/2.0)*(exp_terms[0] + (1.0/De)*( -Gli.Ge*exp_terms[1] - Gri.Ge*exp_terms[2] + 2.0*Gli.Ge*Gri.Ge*cos_term ));
	
	_GIh = (Zi.Yh/2.0)*(exp_terms[0] + (1.0/Dh)*( -Gli.Gh*exp_terms[1] - Gri.Gh*exp_terms[2] + 2.0*Gli.Gh*Gri.Gh*cos_term ));
	
	return;

}


/*! \brief Compute the off-diagonal voltage at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeWVii(double _z, double _zp, std::complex<double> &_WVe, std::complex<double> &_WVh)
{

	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_WVe = (jkzi/2.0/Zi.Ze/De)*( Gli.Ge*exp_terms[1] - Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term );
	
	_WVh = (jkzi/2.0/Zi.Zh/Dh)*( Gli.Gh*exp_terms[1] - Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term );
	
	return;

}


/*! \brief Compute the off-diagonal current at the observation point z (in layer m) due to a unit-strength voltage source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeWIii(double _z, double _zp, std::complex<double> &_WIe, std::complex<double> &_WIh)
{

	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_WIe = (jkzi/2.0/Zi.Ye/De)*( -Gli.Ge*exp_terms[1] + Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term );
	
	_WIh = (jkzi/2.0/Zi.Yh/Dh)*( -Gli.Gh*exp_terms[1] + Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term );

	return;

}


/*! \brief Compute the z- and zp-dependent same-layer exponential terms.*/
void SpectralMGF::ComputeExpTermsGii(double _z, double _zp, std::array<std::complex<double>, 3> &_exp_terms, std::complex<double> &_cos_term, std::complex<double> &_sin_term)
{

	// double zmin = lm->layers[i].zmin;
	// double zmax = lm->layers[i].zmax;
	// double h = lm->layers[i].h;
	double zmin = lm->GetZmin(i);
	double zmax = lm->GetZmax(i);
	double h = lm->GetHeight(i);

	_exp_terms[0] = std::exp(-jkzi*std::abs(_z - _zp));
	_exp_terms[1] = std::exp(-jkzi*(_z + _zp - 2.0*zmin));
	_exp_terms[2] = std::exp(-jkzi*(2.0*zmax - (_z + _zp)));

	_cos_term = (std::exp(jkzi*(-2.0*h + _z - _zp)) + std::exp(-jkzi*(2.0*h + _z - _zp)))/2.0;
	_sin_term = (std::exp(jkzi*(-2.0*h + _z - _zp)) - std::exp(-jkzi*(2.0*h + _z - _zp)))/2.0/J;

	if (i == -1)
	{
		_exp_terms[2] = 0.0;
		_cos_term = 0.0;
		_sin_term = 0.0;
	}
	else if (i == (int)lm->layers.size())
	{
		_exp_terms[1] = 0.0;
		_cos_term = 0.0;
		_sin_term = 0.0;
	}

	exp_terms_computed = true;

	return;

}


/*! \brief Compute the voltage transfer function between the observation point z (in layer m) and the right or left terminal (i.e. upper or lower interface) of layer i, depending on their relative position, for the equivalent transmission line network.*/
void SpectralMGF::ComputeTV(double _z, std::complex<double> &_TVe, std::complex<double> &_TVh)
{

	if (!trf_terms_computed)
		ComputeTrfTerms(_z, trf_terms);

	_TVe = (trf_terms[0] + Gm.Ge*trf_terms[1])*MVe;
	_TVh = (trf_terms[0] + Gm.Gh*trf_terms[1])*MVh;

	return;

}


/*! \brief Compute the current transfer function between the observation point z (in layer m) and the right or left terminal (i.e. upper or lower interface) of layer i, depending on their relative position, for the equivalent transmission line network.*/
void SpectralMGF::ComputeTI(double _z, std::complex<double> &_TIe, std::complex<double> &_TIh)
{

	if (!trf_terms_computed)
		ComputeTrfTerms(_z, trf_terms);

	_TIe = (trf_terms[0] - Gm.Ge*trf_terms[1])*MIe;
	_TIh = (trf_terms[0] - Gm.Gh*trf_terms[1])*MIh;

	return;

}


/*! \brief Compute the z-dependent transfer function terms.*/
void SpectralMGF::ComputeTrfTerms(double _z, std::array<std::complex<double>, 2> &_trf_terms)
{

	// double zmin = lm->layers[m].zmin;
	// double zmax = lm->layers[m].zmax;
	double zmin = lm->GetZmin(m);
	double zmax = lm->GetZmax(m);

	std::fill(_trf_terms.begin(), _trf_terms.end(), 0.0);
	
	if (m < i)
	{
		_trf_terms[0] = std::exp(-jkzm*(_z - zmin));
		_trf_terms[1] = std::exp(-jkzm*(2.0*zmax - _z - zmin));
	}
	else if (m > i)
	{
		_trf_terms[0] = std::exp(-jkzm*(zmax - _z));
		_trf_terms[1] = std::exp(-jkzm*(_z - 2.0*zmin + zmax));
	}
	
	trf_terms_computed = true;

	return;

}


/*! \brief Compute the upwards transfer function multipliers, between the observation layer m and the source layer i of the equivalent transmission line network.*/
void SpectralMGF::ComputeM_Upward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh)
{

	_MVe = 1.0;
	_MVh = 1.0;
	_MIe = 1.0;
	_MIh = 1.0;
	
	for (int ii = m; ii <= i-2; ii++)
	{
		// double hk = lm->layers[ii+1].h;
		double hk = lm->GetHeight(ii+1);
		std::complex<double> kzk = ComputeAxialWaveNumber(lm->k[ii+1], krho);
		std::complex<double> jkzk = J*kzk;

		Gamma Gk;
		ComputeReflectionCoefficient_Upward(krho, ii, Gk);
			     
		_MVe *= (1.0 + Gk.Ge)*std::exp(-jkzk*hk)/(1.0 + Gk.Ge*std::exp(-2.0*jkzk*hk));
		_MVh *= (1.0 + Gk.Gh)*std::exp(-jkzk*hk)/(1.0 + Gk.Gh*std::exp(-2.0*jkzk*hk));

		_MIe *= (1.0 - Gk.Ge)*std::exp(-jkzk*hk)/(1.0 - Gk.Ge*std::exp(-2.0*jkzk*hk));
		_MIh *= (1.0 - Gk.Gh)*std::exp(-jkzk*hk)/(1.0 - Gk.Gh*std::exp(-2.0*jkzk*hk));
	}

	// _MVe *= 1.0/(1.0 + Gm.Ge*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MVh *= 1.0/(1.0 + Gm.Gh*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MIe *= 1.0/(1.0 - Gm.Ge*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MIh *= 1.0/(1.0 - Gm.Gh*std::exp(-jkzm*2.0*lm->layers[m].h));
	_MVe *= 1.0/(1.0 + Gm.Ge*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MVh *= 1.0/(1.0 + Gm.Gh*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MIe *= 1.0/(1.0 - Gm.Ge*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MIh *= 1.0/(1.0 - Gm.Gh*std::exp(-jkzm*2.0*lm->GetHeight(m)));

	return;

}


/*! \brief Compute the downwards transfer function multipliers, between the observation layer m and the source layer i of the equivalent transmission line network.*/
void SpectralMGF::ComputeM_Downward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh)
{

	_MVe = 1.0;
	_MVh = 1.0;
	_MIe = 1.0;
	_MIh = 1.0;
	
	for (int ii = i+1; ii <= m-1; ii++)
	{
		// double hk = lm->layers[ii].h;
		double hk = lm->GetHeight(ii);
		std::complex<double> kzk = ComputeAxialWaveNumber(lm->k[ii], krho);
		std::complex<double> jkzk = J*kzk;

		Gamma Gk;
		ComputeReflectionCoefficient_Downward(krho, ii, Gk);
			     
		_MVe *= (1.0 + Gk.Ge)*std::exp(-jkzk*hk)/(1.0 + Gk.Ge*std::exp(-2.0*jkzk*hk));
		_MVh *= (1.0 + Gk.Gh)*std::exp(-jkzk*hk)/(1.0 + Gk.Gh*std::exp(-2.0*jkzk*hk));

		_MIe *= (1.0 - Gk.Ge)*std::exp(-jkzk*hk)/(1.0 - Gk.Ge*std::exp(-2.0*jkzk*hk));
		_MIh *= (1.0 - Gk.Gh)*std::exp(-jkzk*hk)/(1.0 - Gk.Gh*std::exp(-2.0*jkzk*hk));
	}
	
	// _MVe *= 1.0/(1.0 + Gm.Ge*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MVh *= 1.0/(1.0 + Gm.Gh*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MIe *= 1.0/(1.0 - Gm.Ge*std::exp(-jkzm*2.0*lm->layers[m].h));
	// _MIh *= 1.0/(1.0 - Gm.Gh*std::exp(-jkzm*2.0*lm->layers[m].h));
	_MVe *= 1.0/(1.0 + Gm.Ge*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MVh *= 1.0/(1.0 + Gm.Gh*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MIe *= 1.0/(1.0 - Gm.Ge*std::exp(-jkzm*2.0*lm->GetHeight(m)));
	_MIh *= 1.0/(1.0 - Gm.Gh*std::exp(-jkzm*2.0*lm->GetHeight(m)));

	return;

}


// ==================================================================================
// Computational drivers - curl
// ==================================================================================

/*! \brief Compute the GEM-curl components of the spectral MGF when the source and observation layers are the same.*/
void SpectralMGF::ComputeCurlGEMii(std::array<std::complex<double>, 4> &curlK)
{	

	// ====== Precompute terms ======

	ComputeVVeii(z, zp, VVe);
	ComputeVVhii(z, zp, VVh);
	ComputeVIhii(z, zp, VIh);
	ComputeIVeii(z, zp, IVe);
	
	if (extract_homogeneous && !extract_quasistatic)
	{
		std::complex<double> exp_hom = std::exp(-jkzi*std::abs(z - zp));
		VVe -= (P/2.0)*exp_hom;
		VVh -= (P/2.0)*exp_hom;
		VIh -= (Zi.Zh/2.0)*exp_hom;
		IVe -= (Zi.Ye/2.0)*exp_hom;
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (VVe - VVh)/krhosq;
	}

	if (components_curl[1])
	{
		curlK[1] = VVe + VVh;
	}

	if (components_curl[2])
	{
		curlK[2] = VIh/(J*omega*mui);
	}

	if (components_curl[3])
	{
		curlK[3] = IVe/(J*omega*epsi);
	}
	
	return;

}


/*! \brief Compute the GHJ-curl components of the spectral MGF when the source and observation layers are the same.*/
void SpectralMGF::ComputeCurlGHJii(std::array<std::complex<double>, 4> &curlK)
{	

	// ====== Precompute terms ======

	ComputeIIhii(z, zp, IIh);
	ComputeIIeii(z, zp, IIe);
	ComputeIVeii(z, zp, IVe);
	ComputeVIhii(z, zp, VIh);
	
	if (extract_homogeneous && !extract_quasistatic)
	{
		std::complex<double> exp_hom = std::exp(-jkzi*std::abs(z - zp));
		IIh -= (P/2.0)*exp_hom;
		IIe -= (P/2.0)*exp_hom;
		IVe -= (Zi.Ye/2.0)*exp_hom;
		VIh -= (Zi.Zh/2.0)*exp_hom;
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (IIh - IIe)/krhosq;
	}

	if (components_curl[1])
	{
		curlK[1] = IIh + IIe;
	}

	if (components_curl[2])
	{
		curlK[2] = IVe/(J*omega*epsi);
	}

	if (components_curl[3])
	{
		curlK[3] = VIh/(J*omega*mui);
	}

	return;

}


/*! \brief Compute the GEM-curl of the spectral MGF when the source and observation layers are not the same.*/
void SpectralMGF::ComputeCurlGEMmi(std::array<std::complex<double>, 4> &curlK)
{

	// ====== Precompute terms ======
 
	double _z;
	
	if (m < i)
		// _z = lm->layers[i].zmax;
		_z = lm->GetZmax(i);
	else if (m > i)
		// _z = lm->layers[i].zmin;
		_z = lm->GetZmin(i);

	std::complex<double> GVe_ii, GVh_ii, GIe_ii, GIh_ii;
	ComputeGVii(_z, zp, GVe_ii, GVh_ii);
	ComputeGIii(_z, zp, GIe_ii, GIh_ii);

	ComputeTV(z, TVe, TVh);
	ComputeTI(z, TIe, TIh);

	std::complex<double> GVh_mi = GVh_ii*TVh;
	std::complex<double> GIe_mi = GIe_ii*TIe;

	VIh = GVh_mi;
	IVe = GIe_mi;

	// For VV, the m < i case is treated as the reciprocal to the m > i case. Note that this involves some redundant computations of reflection coefficients etc. and could probably be improved, but has a minimal impact on overall computational cost.

	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver();
		reciprocal = true;
	}

	// _z = lm->layers[i].zmin;
	_z = lm->GetZmin(i);
	
	ComputeTrfTerms(z, trf_terms);

	if (!reciprocal)
	{

		TVe = (trf_terms[0] + Gm.Ge*trf_terms[1])*MVe;
		TVh = (trf_terms[0] + Gm.Gh*trf_terms[1])*MVh;

		std::complex<double> VVeii, VVhii;
		ComputeVVeii(_z, zp, VVeii);
		ComputeVVhii(_z, zp, VVhii);

		VVe = VVeii*TVe;
		VVh = VVhii*TVh;

	}
	else
	{

		TIe = -Zm.Ye*(trf_terms[0] - Gm.Ge*trf_terms[1])*MVe;
		TIh = -Zm.Yh*(trf_terms[0] - Gm.Gh*trf_terms[1])*MVh;

		std::complex<double> VIeii, VIhii;
		ComputeVIeii(_z, zp, VIeii);
		ComputeVIhii(_z, zp, VIhii);

		VVe = -VIeii*TIe;
		VVh = -VIhii*TIh;

		// Reset
		SwapSourceAndObserver();
		reciprocal = false;
		
	}
	
	if (extract_homogeneous && !extract_quasistatic)
	{
		std::complex<double> exp_hom = std::exp(-jkzi*std::abs(z - zp));
		VVe -= (P/2.0)*exp_hom;
		VVh -= (P/2.0)*exp_hom;
		VIh -= (omega*mui/2.0/kzi)*exp_hom;
		IVe -= (omega*epsm/2.0/kzi)*exp_hom;
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (VVe - VVh)/krhosq;
	}

	if (components_curl[1])
	{
		curlK[1] = VVe + VVh;
	}

	if (components_curl[2])
	{
		curlK[2] = VIh/(J*omega*mui);
	}

	if (components_curl[3])
	{
		curlK[3] = IVe/(J*omega*epsm);
	}
	
	return;

}


/*! \brief Compute the GHJ-curl of the spectral MGF when the source and observation layers are not the same.*/
void SpectralMGF::ComputeCurlGHJmi(std::array<std::complex<double>, 4> &curlK)
{

	// ====== Precompute terms ======
 
	double _z;
	
	if (m < i)
		// _z = lm->layers[i].zmax;
		_z = lm->GetZmax(i);
	else if (m > i)
		// _z = lm->layers[i].zmin;
		_z = lm->GetZmin(i);

	std::complex<double> GVe_ii, GVh_ii, GIe_ii, GIh_ii;
	ComputeGVii(_z, zp, GVe_ii, GVh_ii);
	ComputeGIii(_z, zp, GIe_ii, GIh_ii);

	ComputeTV(z, TVe, TVh);
	ComputeTI(z, TIe, TIh);

	std::complex<double> GVh_mi = GVh_ii*TVh;
	std::complex<double> GIe_mi = GIe_ii*TIe;

	VIh = GVh_mi;
	IVe = GIe_mi;

	// For II, the m < i case is treated as the reciprocal to the m > i case. Note that this involves some redundant computations of reflection coefficients etc. and could probably be improved, but has a minimal impact on overall computational cost.

	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver();
		reciprocal = true;
	}

	// _z = lm->layers[i].zmin;
	_z = lm->GetZmin(i);
	
	ComputeTrfTerms(z, trf_terms);

	if (!reciprocal)
	{

		TIe = -Zm.Ye*(trf_terms[0] - Gm.Ge*trf_terms[1])*MVe;
		TIh = -Zm.Yh*(trf_terms[0] - Gm.Gh*trf_terms[1])*MVh;

		std::complex<double> VIeii, VIhii;
		ComputeVIeii(_z, zp, VIeii);
		ComputeVIhii(_z, zp, VIhii);

		VVe = VIeii*TIe;
		VVh = VIhii*TIh;

	}
	else
	{

		TVe = (trf_terms[0] + Gm.Ge*trf_terms[1])*MVe;
		TVh = (trf_terms[0] + Gm.Gh*trf_terms[1])*MVh;

		std::complex<double> VVeii, VVhii;
		ComputeVVeii(_z, zp, VVeii);
		ComputeVVhii(_z, zp, VVhii);

		VVe = -VVeii*TVe;
		VVh = -VVhii*TVh;

		// Reset
		SwapSourceAndObserver();
		reciprocal = false;
		
	}
	
	if (extract_homogeneous && !extract_quasistatic)
	{
		std::complex<double> exp_hom = std::exp(-jkzi*std::abs(z - zp));
		VVe -= (P/2.0)*exp_hom;
		VVh -= (P/2.0)*exp_hom;
		IVe -= (omega*epsi/2.0/kzi)*exp_hom;
		VIh -= (omega*mum/2.0/kzi)*exp_hom;
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (IIh - IIe)/krhosq;
	}

	if (components_curl[1])
	{
		curlK[1] = IIh + IIe;
	}

	if (components_curl[2])
	{
		curlK[2] = IVe/(J*omega*epsi);
	}

	if (components_curl[3])
	{
		curlK[3] = VIh/(J*omega*mum);
	}

	return;

}


// ==================================================================================
// Computational helpers - curl
// ==================================================================================

/*! \brief Compute the TE voltage at the observation point z (in layer m) due to a unit-strength voltage source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeVVhii(double _z, double _zp, std::complex<double> &_VVh)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_VVh = (1.0/2.0)*(P*exp_terms[0] - (1.0/Dh)*( Gli.Gh*exp_terms[1] - Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term ));

	return;

}


/*! \brief Compute the TM voltage at the observation point z (in layer m) due to a unit-strength voltage source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeVVeii(double _z, double _zp, std::complex<double> &_VVe)
{

	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_VVe = (1.0/2.0)*(P*exp_terms[0] - (1.0/De)*( Gli.Ge*exp_terms[1] - Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term ));

	return;

}


/*! \brief Compute the TE current at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeIIhii(double _z, double _zp, std::complex<double> &_IIh)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);
	// ComputeExpTermsGii_dz(_z, _zp, exp_terms, cos_term, sin_term);
		
	_IIh = (1.0/2.0)*(P*exp_terms[0] - (1.0/Dh)*( -Gli.Gh*exp_terms[1] + Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term ));

	return;

}


/*! \brief Compute the TM current at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m.*/
void SpectralMGF::ComputeIIeii(double _z, double _zp, std::complex<double> &_IIe)
{

	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);
	// ComputeExpTermsGii_dz(_z, _zp, exp_terms, cos_term, sin_term);
	
	_IIe = (1.0/2.0)*(P*exp_terms[0] - (1.0/De)*( -Gli.Ge*exp_terms[1] + Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term ));

	return;

}


/*! \brief Compute the TE voltage at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m. Note: this is the same as GVhii, but is coded here separately for readability and modularity.*/
void SpectralMGF::ComputeVIhii(double _z, double _zp, std::complex<double> &_VIh)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_VIh = (Zi.Zh/2.0)*(exp_terms[0] + (1.0/Dh)*( Gli.Gh*exp_terms[1] + Gri.Gh*exp_terms[2] + 2.0*Gli.Gh*Gri.Gh*cos_term ));

	return;

}


/*! \brief Compute the TM voltage at the observation point z (in layer m) due to a unit-strength current source at source point zp (in layer i) of the equivalent transmission line network, for i = m. Note: this is the same as GVeii, but is coded here separately for readability and modularity.*/
void SpectralMGF::ComputeVIeii(double _z, double _zp, std::complex<double> &_VIe)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_VIe = (Zi.Ze/2.0)*(exp_terms[0] + (1.0/De)*( Gli.Ge*exp_terms[1] + Gri.Ge*exp_terms[2] + 2.0*Gli.Ge*Gri.Ge*cos_term ));

	return;

}


/*! \brief Compute the TM current at the observation point z (in layer m) due to a unit-strength voltage source at source point zp (in layer i) of the equivalent transmission line network, for i = m. Note: this is the same as GIeii, but is coded here separately for readability and modularity.*/
void SpectralMGF::ComputeIVeii(double _z, double _zp, std::complex<double> &_IVe)
{
	
	if (!exp_terms_computed)
		ComputeExpTermsGii(_z, _zp, exp_terms, cos_term, sin_term);

	_IVe = (Zi.Ye/2.0)*(exp_terms[0] + (1.0/De)*( -Gli.Ge*exp_terms[1] - Gri.Ge*exp_terms[2] + 2.0*Gli.Ge*Gri.Ge*cos_term ));
	
	return;

}


/*! \brief Compute the z- and zp-dependent same-layer exponential terms in z-derivative form.*/
void SpectralMGF::ComputeExpTermsGii_dz(double _z, double _zp, std::array<std::complex<double>, 3> &_exp_terms, std::complex<double> &_cos_term, std::complex<double> &_sin_term)
{

	// double zmin = lm->layers[i].zmin;
	// double zmax = lm->layers[i].zmax;
	// double h = lm->layers[i].h;
	double zmin = lm->GetZmin(i);
	double zmax = lm->GetZmax(i);
	double h = lm->GetHeight(i);

	if (extract_homogeneous || extract_quasistatic)
		_exp_terms[0] = 0.0;
	else
		_exp_terms[0] = std::exp(-jkzi*std::abs(_z - _zp));
	_exp_terms[1] = std::exp(-jkzi*(_z + _zp - 2.0*zmin))/(-jkzi);
	_exp_terms[2] = std::exp(-jkzi*(2.0*zmax - (_z + _zp)))/(jkzi);

	_cos_term = (std::exp(jkzi*(-2.0*h + _z - _zp))/(jkzi) + std::exp(-jkzi*(2.0*h + _z - _zp))/(-jkzi))/2.0;
	_sin_term = (std::exp(jkzi*(-2.0*h + _z - _zp))/(jkzi) - std::exp(-jkzi*(2.0*h + _z - _zp))/(-jkzi))/2.0/J;

	if (i == -1)
	{
		_exp_terms[2] = 0.0;
		_cos_term = 0.0;
		_sin_term = 0.0;
	}
	else if (i == (int)lm->layers.size())
	{
		_exp_terms[1] = 0.0;
		_cos_term = 0.0;
		_sin_term = 0.0;
	}
	
	exp_terms_computed = false;

	return;

}


/*! \brief Swap the source and observation point and redo all associated precomputations. Useful when using reciprocity to compute components.*/
void SpectralMGF::SwapSourceAndObserver()
{
	
	int i_old = i, m_old = m;
	double z_old = z, zp_old = zp;
	SetLayers(m_old, i_old);
	SetSourcePoint(z_old);
	SetObservationPoint(zp_old);
	SetRadialWaveNumber(krho);
	// P *= -1.0;

	if (z > zp)
		P = 1.0;
	else if (z < zp)
		P = -1.0;
	else
		P = 0.0;	

	return;
	
}



