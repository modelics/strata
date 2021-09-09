/******************************** quasistatic_MGF.cpp ******************************

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


#include <cmath>
#include <complex>
#include <stdexcept>

#include "quasistatic_MGF.hpp"
#include "constants.hpp"

using namespace strata;


// ==================================================================================
// Interface
// ==================================================================================

/*! \brief Initialize the quasistatic MGF engine for a given frequency and layered environment.*/
void QuasistaticMGF::Initialize(LayerManager *_lm, double _f, std::vector<bool> _components, bool _extract_singularities)
{
	
	lm = _lm;
	f = _f;
	omega = 2.0*M_PI*f;

	components = _components;
	extract_singularities = _extract_singularities;
	
	if (!lm->layers_processed)
		lm->ProcessLayers(f);

	// Compute fresnel reflection coefficients for each layer with its adjacent layer
	ComputeFresnelCoefficients();

	initialized = true;
	layers_set = false;
	
	return;
	
}


/*! \brief Set source and observation layers and execute related precomputations.*/
void QuasistaticMGF::SetLayers(int _i, int _m)
{

	// ------ Sanity checks ------
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::SetLayers(): Must call Initialize() before setting layer indices.");
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

	if (m < i)
	{
		// Observation layer is above the source layer
		Rm = GetFresnelCoefficient_Upward(m-1);
		ComputeM_Upward(MVe, MVh, MIe, MIh, nu);
		P = 1.0;
	}
	else if (m > i)
	{
		// Observation layer is below the source layer
		Rm = GetFresnelCoefficient_Downward(m);
		ComputeM_Downward(MVe, MVh, MIe, MIh, nu);
		P = -1.0;
	}
	
	layers_set = true;
	
	return;
	
}


/*! \brief Set the curl of the quasistatic MGF to G_EM, for EFIE-type formulations.*/
void QuasistaticMGF::SetCurlToGEM()
{
	curl_GEM = true;
	return;
}


/*! \brief Set the curl of the quasistatic MGF to G_HJ, for MFIE-type formulations.*/
void QuasistaticMGF::SetCurlToGHJ()
{
	curl_GEM = false;
	return;
}


/*! \brief Driver function to compute desired components of the spectral quasistatic MGF.*/
void QuasistaticMGF::ComputeQMGF_Spectral(double z, double zp, std::complex<double> krho)
{

	// ------ Check that precomputations have been done ------
		
	if (!layers_set)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::ComputeQMGF_Spectral(): Must call SetLayers() before computing the spectral quasistatic MGF.");
	}
	
	// -------------------------------------------------------

	std::fill(K_spectral.begin(), K_spectral.end(), 0.0);
	
	if (i == m)
		ComputeKii_Spectral(z, zp, krho, K_spectral);
	else
		ComputeKmi_Spectral(z, zp, krho, K_spectral);

	K_spectral[0] *= mui/mu0;
	K_spectral[1] *= mui/mu0;
	K_spectral[2] *= mui/mu0;
	K_spectral[3] *= mui/mu0;
	K_spectral[4] *= eps0/epsi;
	
	return;
	
}


/*! \brief Extract computed spectral quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetResult_Spectral(int component)
{
	return K_spectral[component];
}


/*! \brief Driver function to compute desired components of the curl of the spectral quasistatic MGF.*/
void QuasistaticMGF::ComputeCurlQMGF_Spectral(double z, double zp, std::complex<double> krho)
{

	// ------ Check that precomputations have been done ------
		
	if (!layers_set)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::ComputeCurlQMGF_Spectral(): Must call SetLayers() before computing the spectral quasistatic MGF.");
	}
	
	// -------------------------------------------------------

	if (z > zp)
		P = 1.0;
	else if (z < zp)
		P = -1.0;
	else
		P = 0.0;

	std::fill(curlK_spectral.begin(), curlK_spectral.end(), 0.0);
	
	if (curl_GEM)
	{
		if (i == m)
			ComputeCurlGEMii_Spectral(z, zp, krho, curlK_spectral);
		else
			ComputeCurlGEMmi_Spectral(z, zp, krho, curlK_spectral);
	}
	else
	{
		if (i == m)
			ComputeCurlGHJii_Spectral(z, zp, krho, curlK_spectral);
		else
			ComputeCurlGHJmi_Spectral(z, zp, krho, curlK_spectral);
	}
	
	return;
	
}


/*! \brief Extract computed curl of the spectral quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetResultCurl_Spectral(int component)
{
	return curlK_spectral[component];
}


/*! \brief Driver function to compute desired components of the spatial quasistatic MGF.*/
void QuasistaticMGF::ComputeQMGF_Spatial(double z, double zp, double rho)
{

	// ------ Check that precomputations have been done ------

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::ComputeQMGF_Spatial(): Must call SetLayers() before computing the spatial quasistatic MGF.");
	}
	
	// -------------------------------------------------------

	std::fill(K_spatial.begin(), K_spatial.end(), 0.0);
	
	if (i == m)
		ComputeKii_Spatial(z, zp, rho, K_spatial);
	else
		ComputeKmi_Spatial(z, zp, rho, K_spatial);

	K_spatial[0] *= mui/mu0/(4.0*M_PI);
	K_spatial[1] *= mui/mu0/(4.0*M_PI);
	K_spatial[2] *= mui/mu0/(4.0*M_PI);
	K_spatial[3] *= mui/mu0/(4.0*M_PI);
	K_spatial[4] *= eps0/epsi/(4.0*M_PI);

	if (extract_singularities)
		ComputeSingularityFactors(z, zp, rho);

	return;
	
}


/*! \brief Driver function to compute singularity factors.*/
void QuasistaticMGF::ComputeSingularityFactors(double z, double zp, double rho)
{

	// ------ Check that precomputations have been done ------

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::ComputeSingularityFactors(): Must call SetLayers() before computing singularity factors.");
	}
	
	// -------------------------------------------------------

	std::fill(F.begin(), F.end(), 0.0);
	
	if (i == m)
		ComputeFii(z, zp, rho);
	else
		ComputeFmi(z, zp, rho);

	F[0] *= mui/mu0/(4.0*M_PI);
	F[1] *= mui/mu0/(4.0*M_PI);
	F[2] *= mui/mu0/(4.0*M_PI);
	F[3] *= mui/mu0/(4.0*M_PI);
	F[4] *= eps0/epsi/(4.0*M_PI);
	
	return;
	
}


/*! \brief Extract computed spatial quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetResult_Spatial(int component)
{
	return K_spatial[component];
}


/*! \brief Driver function to compute desired components of the curl of the spatial quasistatic MGF.*/
void QuasistaticMGF::ComputeCurlQMGF_Spatial(double z, double zp, double rho)
{

	// ------ Check that precomputations have been done ------
		
	if (!layers_set)
	{
		throw std::logic_error("[ERROR] QuasistaticMGF::ComputeCurlQMGF_Spatial(): Must call SetLayers() before computing the spatial quasistatic MGF.");
	}
	
	// -------------------------------------------------------

	if (z > zp)
		P = 1.0;
	else if (z < zp)
		P = -1.0;
	else
		P = 0.0;
	
	std::fill(curlK_spatial.begin(), curlK_spatial.end(), 0.0);
	std::fill(curlF.begin(), curlF.end(), 0.0);

	if (curl_GEM)
	{
		if (i == m)
			ComputeCurlGEMii_Spatial(z, zp, rho, curlK_spatial);
		else
			ComputeCurlGEMmi_Spatial(z, zp, rho, curlK_spatial);
	}
	else
	{
		if (i == m)
			ComputeCurlGHJii_Spatial(z, zp, rho, curlK_spatial);
		else
			ComputeCurlGHJmi_Spatial(z, zp, rho, curlK_spatial);
	}
	
	curlK_spatial[0] /= 4.0*M_PI;
	curlK_spatial[1] /= 4.0*M_PI;
	curlK_spatial[2] /= 4.0*M_PI;
	curlK_spatial[3] /= 4.0*M_PI;

	return;
	
}


/*! \brief Extract computed curl of the spatial quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetResultCurl_Spatial(int component)
{
	return curlK_spatial[component];
}


/*! \brief Extract computed singularity factor for the quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetSingularityFactor(int component)
{
	return F[component];
}


/*! \brief Extract computed singularity factor for the curl of the quasistatic MGF.*/
std::complex<double> QuasistaticMGF::GetSingularityFactorCurl(int component)
{
	return curlF[component];
}


// ==================================================================================
// Precomputation
// ==================================================================================

/*! \brief Function to compute the upward- and downward-looking Fresnel reflection coefficients for each layer.*/
void QuasistaticMGF::ComputeFresnelCoefficients()
{	
	
	// ====== Upward reflection coefficients ======

	Rrs.clear();
	Rrs.resize(lm->layers.size());

	if (lm->isPEC_top)
	{
		Rrs[0].Fh = -1.0;
		Rrs[0].Fe = -1.0;
	}
	else
	{
		Rrs[0].Fh = -(lm->mu[0] - lm->mu_top)/(lm->mu[0] + lm->mu_top);
		Rrs[0].Fe = (lm->eps[0] - lm->eps_top)/(lm->eps[0] + lm->eps_top);
	}

	for (int ii = 1; ii < (int)lm->layers.size(); ii++)
	{
		Rrs[ii].Fh = -(lm->mu[ii] - lm->mu[ii-1])/(lm->mu[ii] + lm->mu[ii-1]);
		Rrs[ii].Fe = (lm->eps[ii] - lm->eps[ii-1])/(lm->eps[ii] + lm->eps[ii-1]);
	}


	// ====== Downward reflection coefficients ======

	Rls.clear();
	Rls.resize(lm->layers.size());

	if (lm->isPEC_bot)
	{
		Rls.back().Fh = -1.0;
		Rls.back().Fe = -1.0;
	}
	else
	{
		Rls.back().Fh = -(lm->mu.back() - lm->mu_bot)/(lm->mu.back() + lm->mu_bot);
		Rls.back().Fe = (lm->eps.back() - lm->eps_bot)/(lm->eps.back() + lm->eps_bot);
	}

	for (int ii = 0; ii < (int)lm->layers.size()-1; ii++)
	{
		Rls[ii].Fh = -(lm->mu[ii] - lm->mu[ii+1])/(lm->mu[ii] + lm->mu[ii+1]);
		Rls[ii].Fe = (lm->eps[ii] - lm->eps[ii+1])/(lm->eps[ii] + lm->eps[ii+1]);
	}

	return;

}


/*! \brief Function return the upwards Fresnel coefficient at a particular layer interface.*/
Fresnel QuasistaticMGF::GetFresnelCoefficient_Upward(int idx)
{
	if (idx == -2)
	{
		Fresnel Rf;
		Rf.Fh = 0.0;
		Rf.Fe = 0.0;
		return Rf;
	}

	return Rrs[idx+1];
}


/*! \brief Function return the downwards Fresnel coefficient at a particular layer interface.*/
Fresnel QuasistaticMGF::GetFresnelCoefficient_Downward(int idx)
{
	if (idx == (int)(lm->layers.size()))
	{
		Fresnel Rf;
		Rf.Fh = 0.0;
		Rf.Fe = 0.0;
		return Rf;
	}
	
	return Rls[idx];
}


// ==================================================================================
// Helpers
// ==================================================================================

/*! \brief Function to compute kz from k and krho, and make sure it lies in the lower-right complex quadrant.*/
std::complex<double> QuasistaticMGF::ComputeAxialWaveNumber(std::complex<double> k, std::complex<double> krho)
{
	std::complex<double> kz = std::sqrt(k*k - krho*krho);
	CheckQuadrants(kz);
	return kz;
}


/*! \brief Function to ensure that kz lies in the bottom right quadrant of the complex plane.*/
void QuasistaticMGF::CheckQuadrants(std::complex<double> &kz_test)
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

/*! \brief Compute the quasistatic spectral MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeKii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &K)
{

	// ====== Precompute terms ======

	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);
	
	std::array<std::complex<double>, 5> exp_terms_S0, exp_terms_S1;

	if (components[0] || components[3] || components[4])
		ComputeExpTermsKii_kz(z, zp, kzi, exp_terms_S0);
	if (components[1] || components[2])
		ComputeExpTermsKii_krho(z, zp, krho, exp_terms_S1);
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	
	
	// ====== Assemble required component(s) ======

	if (components[0])
	{
		K[0] = (exp_terms_S0[0] + Rl.Fh*exp_terms_S0[1] + Rr.Fh*exp_terms_S0[2] + Rl.Fh*Rr.Fh*(exp_terms_S0[3] + exp_terms_S0[4]))/(2.0*J*kzi);
	}

	if (components[1])
	{
		K[1] = ((Rl.Fe - Rl.Fh)*exp_terms_S1[1] + (-Rr.Fe + Rr.Fh)*exp_terms_S1[2] + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1[3] - exp_terms_S1[4]))/(-2.0*krho*krho);
	}

	if (components[2])
	{
		K[2] = ((-Rl.Fe + Rl.Fh)*exp_terms_S1[1] + (Rr.Fe - Rr.Fh)*exp_terms_S1[2] + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1[3] - exp_terms_S1[4]))/(-2.0*krho*krho);
	}

	if (components[3])
	{
		K[3] = (exp_terms_S0[0] + (-2.0*Rl.Fe + Rl.Fh)*exp_terms_S0[1] + (-2.0*Rr.Fe + Rr.Fh)*exp_terms_S0[2] + (2.0*Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S0[3] + exp_terms_S0[4]))/(2.0*J*kzi);
	}

	if (components[4])
	{
		K[4] = (exp_terms_S0[0] + Rl.Fe*exp_terms_S0[1] + Rr.Fe*exp_terms_S0[2] + Rl.Fe*Rr.Fe*(exp_terms_S0[3] + exp_terms_S0[4]))/(2.0*J*kzi);
	}
		
	return;

}


/*! \brief Compute the quasistatic spectral MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeKmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &K)
{	

	// ====== Precompute terms ======

	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);

	std::array<std::complex<double>, 10> exp_terms_S0, exp_terms_S1;
	ComputeExpTermsKmi_krho(z, zp, krho, exp_terms_S1);

	// If using krho-form
	// exp_terms_S0 = exp_terms_S1;
	// std::complex<double> divisor_S0 = 2.0*krho;

	// If using kz-form
	ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_S0);
	std::complex<double> divisor_S0 = 2.0*J*kzi;

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	
	
	// ====== Assemble required component(s) ======

	if (components[0])
	{
		K[0] = (exp_terms_S0[0] + Rl.Fh*exp_terms_S0[1] + Rr.Fh*exp_terms_S0[2] + Rl.Fh*Rr.Fh*(exp_terms_S0[3] + exp_terms_S0[4]) + Rm.Fh*(exp_terms_S0[5] + Rl.Fh*exp_terms_S0[6] + Rr.Fh*exp_terms_S0[7] + Rl.Fh*Rr.Fh*(exp_terms_S0[8] + exp_terms_S0[9])))*MVh/(divisor_S0);
	}

	if (components[1])
	{
		std::complex<double> c1 = (epsi/epsm)*MVe;
		std::complex<double> c2 = (mum/mui)*MVh;

		K[1] = ((c1 - c2)*exp_terms_S1[0]
				+ (-Rl.Fe*c1 + Rl.Fh*c2)*exp_terms_S1[1]
				+ (-Rr.Fe*c1 + Rr.Fh*c2)*exp_terms_S1[2]
				+ (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[3] + exp_terms_S1[4])
				+ (Rm.Fe*c1 - Rm.Fh*c2)*exp_terms_S1[5]
				+ (-Rm.Fe*Rl.Fe*c1 + Rm.Fh*Rl.Fh*c2)*exp_terms_S1[6]
				+ (-Rm.Fe*Rr.Fe*c1 + Rm.Fh*Rr.Fh*c2)*exp_terms_S1[7]
				+ (Rm.Fe*Rl.Fe*Rr.Fe*c1 - Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[8] + exp_terms_S1[9]))*P/(2.0*krho*krho);
	}

	if (components[2])
	{
		std::complex<double> c1 = (epsm*mum/epsi/mui)*MIe;
		std::complex<double> c2 = MIh;

		K[2] = ((c1 - c2)*exp_terms_S1[0]
				+ (Rl.Fe*c1 - Rl.Fh*c2)*exp_terms_S1[1]
				+ (Rr.Fe*c1 - Rr.Fh*c2)*exp_terms_S1[2]
				+ (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[3] + exp_terms_S1[4])
				+ (-Rm.Fe*c1 + Rm.Fh*c2)*exp_terms_S1[5]
				+ (-Rm.Fe*Rl.Fe*c1 + Rm.Fh*Rl.Fh*c2)*exp_terms_S1[6]
				+ (-Rm.Fe*Rr.Fe*c1 + Rm.Fh*Rr.Fh*c2)*exp_terms_S1[7]
				+ (-Rm.Fe*Rl.Fe*Rr.Fe*c1 + Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[8] + exp_terms_S1[9]))*P/(2.0*krho*krho);
	}

	if (components[3])
	{
		// std::complex<double> c = (1.0 + std::pow(ki, 2)/std::pow(km, 2))*MIe;
		std::complex<double> c = (1.0 + epsi*mui/epsm/mum)*MIe;

		K[3] = (((c - MIh)*exp_terms_S0[0]
				 + (-Rl.Fe*c + Rl.Fh*MIh)*exp_terms_S0[1]
				 + (-Rr.Fe*c + Rr.Fh*MIh)*exp_terms_S0[2]
				 + (Rl.Fe*Rr.Fe*c - Rl.Fh*Rr.Fh*MIh)*(exp_terms_S0[3] + exp_terms_S0[4])
				 + (-Rm.Fe*c + Rm.Fh*MIh)*exp_terms_S0[5]
				 + (Rm.Fe*Rl.Fe*c - Rm.Fh*Rl.Fh*MIh)*exp_terms_S0[6]
				 + (Rm.Fe*Rr.Fe*c - Rm.Fh*Rr.Fh*MIh)*exp_terms_S0[7]
				 + (-Rm.Fe*Rl.Fe*Rr.Fe*c + Rm.Fh*Rl.Fh*Rr.Fh*MIh)*(exp_terms_S0[8] + exp_terms_S0[9]))/(divisor_S0))*(mum/mui);
	}

	if (components[4])
	{
		K[4] = (exp_terms_S0[0] + Rl.Fe*exp_terms_S0[1] + Rr.Fe*exp_terms_S0[2] + Rl.Fe*Rr.Fe*(exp_terms_S0[3] + exp_terms_S0[4]) + Rm.Fe*(exp_terms_S0[5] + Rl.Fe*exp_terms_S0[6] + Rr.Fe*exp_terms_S0[7] + Rl.Fe*Rr.Fe*(exp_terms_S0[8] + exp_terms_S0[9])))*MVe/(divisor_S0);
	}
	
	return;

}


/*! \brief Compute the quasistatic spatial MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeKii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 5> &K)
{
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);


	// ====== Compute singularity factors ======

	if (extract_singularities)
	{
		ComputeFii(z, zp, rho);

		// if (std::abs(ki) < 1.0e-15)
		// {
		// 	std::fill(K.begin(), K.end(), 0.0);
		// 	return;
		// }
	}
	

	// ====== Precompute terms ======

	std::array<std::complex<double>, 5> exp_terms_S0, exp_terms_S1;

	if (components[0] || components[3] || components[4])
		ComputeExpTermsKii_SpatialS0(z, zp, rho, exp_terms_S0);
	if (components[1] || components[2])
		ComputeExpTermsKii_SpatialS1(z, zp, rho, exp_terms_S1);

	
	// ====== Assemble required component(s) ======

	if (components[0])
	{
		K[0] = exp_terms_S0[0] + Rl.Fh*exp_terms_S0[1] + Rr.Fh*exp_terms_S0[2] + Rl.Fh*Rr.Fh*(exp_terms_S0[3] + exp_terms_S0[4]);
	}

	if (components[1])
	{
		K[1] = -1.0*((Rl.Fe - Rl.Fh)*exp_terms_S1[1] + (-Rr.Fe + Rr.Fh)*exp_terms_S1[2] + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1[3] - exp_terms_S1[4]));
	}

	if (components[2])
	{
		K[2] = -1.0*((-Rl.Fe + Rl.Fh)*exp_terms_S1[1] + (Rr.Fe - Rr.Fh)*exp_terms_S1[2] + (Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S1[3] - exp_terms_S1[4]));
	}

	if (components[3])
	{
		K[3] = exp_terms_S0[0] + (-2.0*Rl.Fe + Rl.Fh)*exp_terms_S0[1] + (-2.0*Rr.Fe + Rr.Fh)*exp_terms_S0[2] + (2.0*Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_S0[3] + exp_terms_S0[4]);
	}

	if (components[4])
	{
		K[4] = exp_terms_S0[0] + Rl.Fe*exp_terms_S0[1] + Rr.Fe*exp_terms_S0[2] + Rl.Fe*Rr.Fe*(exp_terms_S0[3] + exp_terms_S0[4]);
	}
		
	return;

}


/*! \brief Compute the quasistatic spatial MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeKmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 5> &K)
{	
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);


	// ====== Compute singularity factors ======

	if (extract_singularities)
	{
		ComputeFmi(z, zp, rho);

		// if (std::abs(ki) < 1.0e-15)
		// {
		// 	std::fill(K.begin(), K.end(), 0.0);
		// 	return;
		// }
	}


	// ====== Precompute terms ======

	std::array<std::complex<double>, 10> exp_terms_S0, exp_terms_S1;

	if (components[0] || components[3] || components[4])
		ComputeExpTermsKmi_SpatialS0(z, zp, rho, exp_terms_S0);
	if (components[1] || components[2])
		ComputeExpTermsKmi_SpatialS1(z, zp, rho, exp_terms_S1);

	
	// ====== Assemble required component(s) ======

	if (components[0])
	{
		K[0] = (exp_terms_S0[0] + Rl.Fh*exp_terms_S0[1] + Rr.Fh*exp_terms_S0[2] + Rl.Fh*Rr.Fh*(exp_terms_S0[3] + exp_terms_S0[4]) + Rm.Fh*(exp_terms_S0[5] + Rl.Fh*exp_terms_S0[6] + Rr.Fh*exp_terms_S0[7] + Rl.Fh*Rr.Fh*(exp_terms_S0[8] + exp_terms_S0[9])))*MVh;
	}

	if (components[1])
	{
		std::complex<double> c1 = (epsi/epsm)*MVe;
		std::complex<double> c2 = (mum/mui)*MVh;

		K[1] = ((c1 - c2)*exp_terms_S1[0]
			+ (-Rl.Fe*c1 + Rl.Fh*c2)*exp_terms_S1[1]
			+ (-Rr.Fe*c1 + Rr.Fh*c2)*exp_terms_S1[2]
			+ (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[3] + exp_terms_S1[4])
			+ (Rm.Fe*c1 - Rm.Fh*c2)*exp_terms_S1[5]
			+ (-Rm.Fe*Rl.Fe*c1 + Rm.Fh*Rl.Fh*c2)*exp_terms_S1[6]
			+ (-Rm.Fe*Rr.Fe*c1 + Rm.Fh*Rr.Fh*c2)*exp_terms_S1[7]
				+ (Rm.Fe*Rl.Fe*Rr.Fe*c1 - Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[8] + exp_terms_S1[9]))*P;
	}

	if (components[2])
	{
		std::complex<double> c1 = (epsm*mum/epsi/mui)*MIe;
		std::complex<double> c2 = MIh;

		K[2] = ((c1 - c2)*exp_terms_S1[0]
			+ (Rl.Fe*c1 - Rl.Fh*c2)*exp_terms_S1[1]
			+ (Rr.Fe*c1 - Rr.Fh*c2)*exp_terms_S1[2]
			+ (Rl.Fe*Rr.Fe*c1 - Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[3] + exp_terms_S1[4])
			+ (-Rm.Fe*c1 + Rm.Fh*c2)*exp_terms_S1[5]
			+ (-Rm.Fe*Rl.Fe*c1 + Rm.Fh*Rl.Fh*c2)*exp_terms_S1[6]
			+ (-Rm.Fe*Rr.Fe*c1 + Rm.Fh*Rr.Fh*c2)*exp_terms_S1[7]
				+ (-Rm.Fe*Rl.Fe*Rr.Fe*c1 + Rm.Fh*Rl.Fh*Rr.Fh*c2)*(exp_terms_S1[8] + exp_terms_S1[9]))*P;
	}

	if (components[3])
	{
		// std::complex<double> c = (1.0 + std::pow(ki, 2)/std::pow(km, 2))*MIe;
		std::complex<double> c = (1.0 + epsi*mui/epsm/mum)*MIe;

		K[3] = ((c - MIh)*exp_terms_S0[0]
				+ (-Rl.Fe*c + Rl.Fh*MIh)*exp_terms_S0[1]
				+ (-Rr.Fe*c + Rr.Fh*MIh)*exp_terms_S0[2]
				+ (Rl.Fe*Rr.Fe*c - Rl.Fh*Rr.Fh*MIh)*(exp_terms_S0[3] + exp_terms_S0[4])
				+ (-Rm.Fe*c + Rm.Fh*MIh)*exp_terms_S0[5]
				+ (Rm.Fe*Rl.Fe*c - Rm.Fh*Rl.Fh*MIh)*exp_terms_S0[6]
				+ (Rm.Fe*Rr.Fe*c - Rm.Fh*Rr.Fh*MIh)*exp_terms_S0[7]
				+ (-Rm.Fe*Rl.Fe*Rr.Fe*c + Rm.Fh*Rl.Fh*Rr.Fh*MIh)*(exp_terms_S0[8] + exp_terms_S0[9]))*(mum/mui);
	}

	if (components[4])
	{
		K[4] = (exp_terms_S0[0] + Rl.Fe*exp_terms_S0[1] + Rr.Fe*exp_terms_S0[2] + Rl.Fe*Rr.Fe*(exp_terms_S0[3] + exp_terms_S0[4]) + Rm.Fe*(exp_terms_S0[5] + Rl.Fe*exp_terms_S0[6] + Rr.Fe*exp_terms_S0[7] + Rl.Fe*Rr.Fe*(exp_terms_S0[8] + exp_terms_S0[9])))*MVe;
	}
	
	return;

}


/*! \brief Compute the singularity factors when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeFii(double z, double zp, double rho)
{
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);

	F[0] = 1.0;
	F[3] = 1.0;
	F[4] = 1.0;

	F[1] = 0.0;
	F[2] = 0.0;

	// Source and observer are at the bottom interface of the layer
	// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
	if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
	{
		F[0] += Rl.Fh;
		F[3] += -2.0*Rl.Fe + Rl.Fh;
		F[4] += Rl.Fe;

		F[1] = -1.0*(Rl.Fe - Rl.Fh);
		F[2] = -1.0*(-Rl.Fe + Rl.Fh);
	}
	// Source and observer are at the top interface of the layer
	// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
	else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
	{
		F[0] += Rr.Fh;
		F[3] += -2.0*Rr.Fe + Rr.Fh;
		F[4] += Rr.Fe;

		F[1] = -1.0*(-Rr.Fe + Rr.Fh);
		F[2] = -1.0*(Rr.Fe - Rr.Fh);
	}

	// If rho ~ 0, J_1(rho*krho) ~ 0, so off-diagonal terms are ~0.
	if (LaterallyCoincident(rho))
	{
		F[1] = 0.0;
		F[2] = 0.0;
	}

	if (!components[0])
		F[0] = 0.0;
	if (!components[1])
		F[1] = 0.0;
	if (!components[2])
		F[2] = 0.0;
	if (!components[3])
		F[3] = 0.0;
	if (!components[4])
		F[4] = 0.0;

	return;

}


/*! \brief Compute the singularity factors when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeFmi(double z, double zp, double rho)
{	
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);

	F[0] = MVh;
	// F[3] = (1.0 + std::pow(ki, 2)/std::pow(km, 2))*MIe - MIh;
	F[3] = (1.0 + epsi*mui/epsm/mum)*MIe - MIh;
	F[4] = MVe;

	F[1] = 0.0;
	F[2] = 0.0;

	if (AxiallyCoincident(z, zp))
	{
		F[1] = ((epsi/epsm)*MVe - (mum/mui)*MVh)*P;
		F[2] = ((epsm*mum/epsi/mui)*MIe - MIh)*P;
			
		if (m == i + 1)
		{
			F[0] = (1.0 + Rl.Fh)*MVh;
			// F[3] = (1.0 - Rl.Fe)*(1.0 + std::pow(ki, 2)/std::pow(km, 2))*MIe + (-1.0 + Rl.Fh)*MIh;
			F[3] = (1.0 - Rl.Fe)*(1.0 + epsi*mui/epsm/mum)*MIe + (-1.0 + Rl.Fh)*MIh;
			F[4] = (1.0 + Rl.Fe)*MVe;

			F[1] += (-Rl.Fe*(epsi/epsm)*MVe + Rl.Fh*(mum/mui)*MVh)*P;
			F[2] += (Rl.Fe*(epsm*mum/epsi/mui)*MIe - Rl.Fh*MIh)*P;
		}
		else if (m == i - 1)
		{
			F[0] = (1.0 + Rr.Fh)*MVh;
			// F[3] = (1.0 - Rr.Fe)*(1.0 + std::pow(ki, 2)/std::pow(km, 2))*MIe + (-1.0 + Rr.Fh)*MIh;
			F[3] = (1.0 - Rr.Fe)*(1.0 + epsi*mui/epsm/mum)*MIe + (-1.0 + Rr.Fh)*MIh;
			F[4] = (1.0 + Rr.Fe)*MVe;

			F[1] += (-Rr.Fe*(epsi/epsm)*MVe + Rr.Fh*(mum/mui)*MVh)*P;
			F[2] += (Rr.Fe*(epsm*mum/epsi/mui)*MIe - Rr.Fh*MIh)*P;
		}
	}

	F[3] *= (mum/mui);

	// If rho ~ 0, J_1(rho*krho) ~ 0, so off-diagonal terms are ~0.
	if (LaterallyCoincident(rho))
	{
		F[1] = 0.0;
		F[2] = 0.0;
	}

	if (!components[0])
		F[0] = 0.0;
	if (!components[1])
		F[1] = 0.0;
	if (!components[2])
		F[2] = 0.0;
	if (!components[3])
		F[3] = 0.0;
	if (!components[4])
		F[4] = 0.0;

	return;

}


// ==================================================================================
// Computational drivers - curl
// ==================================================================================

/*! \brief Compute the GEM-curl of the quasistatic spectral MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeCurlGEMii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK)
{

	// ====== Precompute terms ======

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);

	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);
	std::array<std::complex<double>, 5> exp_terms_krho, exp_terms_kz;

	if (components_curl[0])
		ComputeExpTermsKii_krho(z, zp, krho, exp_terms_krho);
	if (components_curl[1] || components_curl[2] || components_curl[3])
		ComputeExpTermsKii_kz(z, zp, kzi, exp_terms_kz);

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = ( (-Rl.Fe + Rl.Fh)*exp_terms_krho[1] + (Rr.Fe - Rr.Fh)*exp_terms_krho[2] + (-Rl.Fe*Rr.Fe + Rl.Fh*Rr.Fh)*(exp_terms_krho[3] - exp_terms_krho[4]) )/(2.0*krho*krho);
	}

	if (components_curl[1])
	{
		curlK[1] = P*exp_terms_kz[0] + ( (-Rl.Fe - Rl.Fh)*exp_terms_kz[1] + (Rr.Fe + Rr.Fh)*exp_terms_kz[2] + (-Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(exp_terms_kz[3] - exp_terms_kz[4]) )/(2.0);
	}

	if (components_curl[2])
	{
		curlK[2] = ( exp_terms_kz[0] + Rl.Fh*exp_terms_kz[1] + Rr.Fh*exp_terms_kz[2] + Rl.Fh*Rr.Fh*(exp_terms_kz[3] + exp_terms_kz[4]) )/(2.0*J*kzi);
	}

	if (components_curl[3])
	{
		curlK[3] = ( exp_terms_kz[0] - Rl.Fe*exp_terms_kz[1] - Rr.Fe*exp_terms_kz[2] + Rl.Fe*Rr.Fe*(exp_terms_kz[3] + exp_terms_kz[4]) )/(2.0*J*kzi);
	}
		
	return;

}


/*! \brief Compute the GHJ-curl of the quasistatic spectral MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeCurlGHJii_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK)
{

	// ====== Precompute terms ======

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);

	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);
	std::array<std::complex<double>, 5> exp_terms_krho, exp_terms_kz;

	if (components_curl[0])
		ComputeExpTermsKii_krho(z, zp, krho, exp_terms_krho);
	if (components_curl[1] || components_curl[2] || components_curl[3])
		ComputeExpTermsKii_kz(z, zp, kzi, exp_terms_kz);


	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = ( (Rl.Fh - Rl.Fe)*exp_terms_krho[1] + (-Rr.Fh + Rr.Fe)*exp_terms_krho[2] + (-Rl.Fh*Rr.Fh + Rl.Fe*Rr.Fe)*(exp_terms_krho[3] - exp_terms_krho[4]) )/(2.0*krho*krho);
	}

	if (components_curl[1])
	{
		curlK[1] = P*exp_terms_kz[0] + ( (Rl.Fh + Rl.Fe)*exp_terms_kz[1] + (-Rr.Fh - Rr.Fe)*exp_terms_kz[2] + (-Rl.Fh*Rr.Fh - Rl.Fe*Rr.Fe)*(exp_terms_kz[3] - exp_terms_kz[4]) )/(2.0);
	}

	if (components_curl[2])
	{
		curlK[2] = ( exp_terms_kz[0] - Rl.Fe*exp_terms_kz[1] - Rr.Fe*exp_terms_kz[2] + Rl.Fe*Rr.Fe*(exp_terms_kz[3] + exp_terms_kz[4]) )/(2.0*J*kzi);
	}
	
	if (components_curl[3])
	{
		curlK[3] = ( exp_terms_kz[0] + Rl.Fh*exp_terms_kz[1] + Rr.Fh*exp_terms_kz[2] + Rl.Fh*Rr.Fh*(exp_terms_kz[3] + exp_terms_kz[4]) )/(2.0*J*kzi);
	}
	
	return;

}


/*! \brief Compute the GEM-curl of the quasistatic spectral MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeCurlGEMmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK)
{	

	// ====== Precompute terms ======

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);
	
	std::array<std::complex<double>, 10> exp_terms_krho, exp_terms_kz;
	std::complex<double> VIh, IVe;

	if (components_curl[2] || components_curl[3])
	{
		
		ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
		
		VIh = (Rl.Fh*exp_terms_kz[1] + Rr.Fh*exp_terms_kz[2] + Rl.Fh*Rr.Fh*(exp_terms_kz[3] + exp_terms_kz[4]) + Rm.Fh*(exp_terms_kz[5] + Rl.Fh*exp_terms_kz[6] + Rr.Fh*exp_terms_kz[7] + Rl.Fh*Rr.Fh*(exp_terms_kz[8] + exp_terms_kz[9])))*MVh/(2.0*J*kzi);

		IVe = ((epsi/epsm)*(-Rl.Fe*exp_terms_kz[1] - Rr.Fe*exp_terms_kz[2] + Rl.Fe*Rr.Fe*(exp_terms_kz[3] + exp_terms_kz[4]) - Rm.Fe*(exp_terms_kz[5] - Rl.Fe*exp_terms_kz[6] - Rr.Fe*exp_terms_kz[7] + Rl.Fe*Rr.Fe*(exp_terms_kz[8] + exp_terms_kz[9]))))*MIe/(2.0*J*kzi);

		std::complex<double> exp_hom = std::exp(-J*kzi*std::abs(z - zp));
		VIh += exp_hom/(2.0*J*kzi);
		IVe += exp_hom/(2.0*J*kzi);

	}
	
	// For VV, the m < i case is treated as the reciprocal to the m > i case.
	
	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
	    kzi = ComputeAxialWaveNumber(ki, krho);
		reciprocal = true;
	}
	
	std::complex<double> VVe_S2, VVe_S0, VVh_S2, VVh_S0;
	
	if (!reciprocal)
	{

		auto ComputeVVe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fe*terms[1] + Rr.Fe*terms[2] - Rl.Fe*Rr.Fe*(terms[3] - terms[4]) + this->Rm.Fe*(this->P*terms[5] - Rl.Fe*terms[6] + Rr.Fe*terms[7] - Rl.Fe*Rr.Fe*(terms[8] - terms[9])))*MVe;
			};

		auto ComputeVVh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fh*terms[1] + Rr.Fh*terms[2] - Rl.Fh*Rr.Fh*(terms[3] - terms[4]) + this->Rm.Fh*(this->P*terms[5] - Rl.Fh*terms[6] + Rr.Fh*terms[7] - Rl.Fh*Rr.Fh*(terms[8] - terms[9])))*MVh;
			};

		if (components_curl[0])
		{
			ComputeExpTermsKmi_krho(z, zp, krho, exp_terms_krho);
			VVe_S2 = ComputeVVe(exp_terms_krho);
			VVh_S2 = ComputeVVh(exp_terms_krho);
		}

		if (components_curl[1])
		{
			ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
			VVe_S0 = ComputeVVe(exp_terms_kz);
			VVh_S0 = ComputeVVh(exp_terms_kz);
		}
		
	}
	else
	{

		auto ComputeIIe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((epsm/epsi)*(Rl.Fe*terms[1] + Rr.Fe*terms[2] + Rl.Fe*Rr.Fe*(terms[3] + terms[4]) - this->Rm.Fe*(terms[5] + Rl.Fe*terms[6] + Rr.Fe*terms[7] + Rl.Fe*Rr.Fe*(terms[8] + terms[9]))))*MVe;
			};

		auto ComputeIIh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((mui/mum)*(Rl.Fh*terms[1] + Rr.Fh*terms[2] + Rl.Fh*Rr.Fh*(terms[3] + terms[4]) - this->Rm.Fh*(terms[5] + Rl.Fh*terms[6] + Rr.Fh*terms[7] + Rl.Fh*Rr.Fh*(terms[8] + terms[9]))))*MVh;
			};

		if (components_curl[0])
		{
			ComputeExpTermsKmi_krho(z, zp, krho, exp_terms_krho);
			VVe_S2 = -ComputeIIe(exp_terms_krho);
			VVh_S2 = -ComputeIIh(exp_terms_krho);
		}

		if (components_curl[1])
		{
			ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
			VVe_S0 = -ComputeIIe(exp_terms_kz);
			VVh_S0 = -ComputeIIh(exp_terms_kz);
		}
		
		// Reset
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		kzi = ComputeAxialWaveNumber(ki, krho);
		reciprocal = false;
		
	}

	std::complex<double> exp_hom_S0 = std::exp(-J*kzi*std::abs(z - zp));
	VVe_S0 += (P/2.0)*exp_hom_S0;
	VVh_S0 += (P/2.0)*exp_hom_S0;


	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (VVe_S2 - VVh_S2)/krho/krho;
	}

	if (components_curl[1])
	{
		curlK[1] = VVe_S0 + VVh_S0;
	}

	if (components_curl[2])
	{
		curlK[2] = VIh;
	}

	if (components_curl[3])
	{
		curlK[3] = IVe;
	}

	return;

}


/*! \brief Compute the GHJ-curl of the quasistatic spectral MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeCurlGHJmi_Spectral(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 4> &curlK)
{	

	// ====== Precompute terms ======

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	std::complex<double> kzi = ComputeAxialWaveNumber(ki, krho);
	
	std::array<std::complex<double>, 10> exp_terms_krho, exp_terms_kz;
	std::complex<double> VIh, IVe;

	if (components_curl[2] || components_curl[3])
	{
		
		ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
		
		VIh = (Rl.Fh*exp_terms_kz[1] + Rr.Fh*exp_terms_kz[2] + Rl.Fh*Rr.Fh*(exp_terms_kz[3] + exp_terms_kz[4]) + Rm.Fh*(exp_terms_kz[5] + Rl.Fh*exp_terms_kz[6] + Rr.Fh*exp_terms_kz[7] + Rl.Fh*Rr.Fh*(exp_terms_kz[8] + exp_terms_kz[9])))*MVh/(2.0*J*kzi);

		IVe = ((epsi/epsm)*(-Rl.Fe*exp_terms_kz[1] - Rr.Fe*exp_terms_kz[2] + Rl.Fe*Rr.Fe*(exp_terms_kz[3] + exp_terms_kz[4]) - Rm.Fe*(exp_terms_kz[5] - Rl.Fe*exp_terms_kz[6] - Rr.Fe*exp_terms_kz[7] + Rl.Fe*Rr.Fe*(exp_terms_kz[8] + exp_terms_kz[9]))))*MIe/(2.0*J*kzi);

		std::complex<double> exp_hom = std::exp(-J*kzi*std::abs(z - zp));
		VIh += exp_hom/(2.0*J*kzi);
		IVe += exp_hom/(2.0*J*kzi);

	}
	
	// For VV, the m < i case is treated as the reciprocal to the m > i case.
	
	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
	    kzi = ComputeAxialWaveNumber(ki, krho);
		reciprocal = true;
	}
	
	std::complex<double> IIe_S2, IIe_S0, IIh_S2, IIh_S0;
	
	if (!reciprocal)
	{

		auto ComputeIIe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((epsm/epsi)*(Rl.Fe*terms[1] + Rr.Fe*terms[2] + Rl.Fe*Rr.Fe*(terms[3] + terms[4]) - this->Rm.Fe*(terms[5] + Rl.Fe*terms[6] + Rr.Fe*terms[7] + Rl.Fe*Rr.Fe*(terms[8] + terms[9]))))*MVe;
			};

		auto ComputeIIh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((mui/mum)*(Rl.Fh*terms[1] + Rr.Fh*terms[2] + Rl.Fh*Rr.Fh*(terms[3] + terms[4]) - this->Rm.Fh*(terms[5] + Rl.Fh*terms[6] + Rr.Fh*terms[7] + Rl.Fh*Rr.Fh*(terms[8] + terms[9]))))*MVh;
			};

		if (components_curl[0])
		{
			ComputeExpTermsKmi_krho(z, zp, krho, exp_terms_krho);
			IIe_S2 = ComputeIIe(exp_terms_krho);
			IIh_S2 = ComputeIIh(exp_terms_krho);
		}

		if (components_curl[1])
		{
			ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
			IIe_S0 = ComputeIIe(exp_terms_kz);
			IIh_S0 = ComputeIIh(exp_terms_kz);
		}
		
	}
	else
	{

		auto ComputeVVe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fe*terms[1] + Rr.Fe*terms[2] - Rl.Fe*Rr.Fe*(terms[3] - terms[4]) + this->Rm.Fe*(this->P*terms[5] - Rl.Fe*terms[6] + Rr.Fe*terms[7] - Rl.Fe*Rr.Fe*(terms[8] - terms[9])))*MVe;
			};

		auto ComputeVVh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fh*terms[1] + Rr.Fh*terms[2] - Rl.Fh*Rr.Fh*(terms[3] - terms[4]) + this->Rm.Fh*(this->P*terms[5] - Rl.Fh*terms[6] + Rr.Fh*terms[7] - Rl.Fh*Rr.Fh*(terms[8] - terms[9])))*MVh;
			};

		if (components_curl[0])
		{
			ComputeExpTermsKmi_krho(z, zp, krho, exp_terms_krho);
			IIe_S2 = -ComputeVVe(exp_terms_krho);
			IIh_S2 = -ComputeVVh(exp_terms_krho);
		}

		if (components_curl[1])
		{
			ComputeExpTermsKmi_kz(z, zp, kzi, exp_terms_kz);
			IIe_S0 = -ComputeVVe(exp_terms_kz);
			IIh_S0 = -ComputeVVh(exp_terms_kz);
		}
		
		// Reset
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		kzi = ComputeAxialWaveNumber(ki, krho);
		reciprocal = false;
		
	}

	std::complex<double> exp_hom_S0 = std::exp(-J*kzi*std::abs(z - zp));
	IIe_S0 += (P/2.0)*exp_hom_S0;
	IIh_S0 += (P/2.0)*exp_hom_S0;


	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (IIh_S2 - IIe_S2)/krho/krho;
	}

	if (components_curl[1])
	{
		curlK[1] = IIh_S0 + IIe_S0;
	}

	if (components_curl[2])
	{
		curlK[2] = IVe;
	}

	if (components_curl[3])
	{
		curlK[3] = VIh;
	}

	return;

}


/*! \brief Compute the GEM-curl of the quasistatic spatial MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeCurlGEMii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK)
{

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	

	// ====== Compute singularity factors ======

	// Note: it is assumed that the homogeneous terms is always extracted
	if (extract_singularities)
	{

		std::fill(curlF.begin(), curlF.end(), 0.0);

		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			curlF[1] = (-Rl.Fe - Rl.Fh)/2.0;
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			curlF[1] = (Rr.Fe + Rr.Fh)/2.0;
		}

	}


	// ====== Precompute terms ======

	std::array<std::complex<double>, 5> terms_S0, terms_S1, terms_S2;

	if (components_curl[0])
		ComputeTermsCurlKii_SpatialS2(z, zp, rho, terms_S2);
	if (components_curl[1])
		ComputeTermsCurlKii_SpatialS0(z, zp, rho, terms_S0);
	if (components_curl[2] || components_curl[3])
		ComputeTermsCurlKii_SpatialS1(z, zp, rho, terms_S1);

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (-Rl.Fe + Rl.Fh)*terms_S2[1] + (Rr.Fe - Rr.Fh)*terms_S2[2] + (-Rl.Fe*Rr.Fe + Rl.Fh*Rr.Fh)*(terms_S2[3] - terms_S2[4]);
	}

	if (components_curl[1])
	{
		curlK[1] = (-Rl.Fe - Rl.Fh)*terms_S0[1] + (Rr.Fe + Rr.Fh)*terms_S0[2] + (-Rl.Fe*Rr.Fe - Rl.Fh*Rr.Fh)*(terms_S0[3] - terms_S0[4]);
	}

	if (components_curl[2])
	{
		curlK[2] = Rl.Fh*terms_S1[1] + Rr.Fh*terms_S1[2] + Rl.Fh*Rr.Fh*(terms_S1[3] + terms_S1[4]);
	}

	if (components_curl[3])
	{
		curlK[3] = -Rl.Fe*terms_S1[1] - Rr.Fe*terms_S1[2] + Rl.Fe*Rr.Fe*(terms_S1[3] + terms_S1[4]);
	}
	
	return;

}


/*! \brief Compute the GHJ-curl of the quasistatic spatial MGF when the source and observation layers are the same.*/
void QuasistaticMGF::ComputeCurlGHJii_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK)
{

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);
	

	// ====== Compute singularity factors ======

	// Note: it is assumed that the homogeneous terms is always extracted
	if (extract_singularities)
	{

		std::fill(curlF.begin(), curlF.end(), 0.0);

		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			curlF[1] = (Rl.Fe + Rl.Fh)/2.0;
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			curlF[1] = (-Rr.Fe - Rr.Fh)/2.0;
		}

	}


	// ====== Precompute terms ======

	std::array<std::complex<double>, 5> terms_S0, terms_S1, terms_S2;

	if (components_curl[0])
		ComputeTermsCurlKii_SpatialS2(z, zp, rho, terms_S2);
	if (components_curl[1])
		ComputeTermsCurlKii_SpatialS0(z, zp, rho, terms_S0);
	if (components_curl[2] || components_curl[3])
		ComputeTermsCurlKii_SpatialS1(z, zp, rho, terms_S1);

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = (Rl.Fh - Rl.Fe)*terms_S2[1] + (-Rr.Fh + Rr.Fe)*terms_S2[2] + (-Rl.Fh*Rr.Fh + Rl.Fe*Rr.Fe)*(terms_S2[3] - terms_S2[4]);
	}

	if (components_curl[1])
	{
		curlK[1] = (Rl.Fh + Rl.Fe)*terms_S0[1] + (-Rr.Fh - Rr.Fe)*terms_S0[2] + (-Rl.Fh*Rr.Fh - Rl.Fe*Rr.Fe)*(terms_S0[3] - terms_S0[4]);
	}

	if (components_curl[2])
	{
		curlK[2] = -Rl.Fe*terms_S1[1] - Rr.Fe*terms_S1[2] + Rl.Fe*Rr.Fe*(terms_S1[3] + terms_S1[4]);
	}

	if (components_curl[3])
	{
		curlK[3] = Rl.Fh*terms_S1[1] + Rr.Fh*terms_S1[2] + Rl.Fh*Rr.Fh*(terms_S1[3] + terms_S1[4]);
	}
	
	return;

}


/*! \brief Compute the GEM-curl of the quasistatic spatial MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeCurlGEMmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK)
{

	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);

	// ====== Compute singularity factors ======

	// Note: it is assumed that the homogeneous terms is always extracted
	if (extract_singularities)
	{

		std::fill(curlF.begin(), curlF.end(), 0.0);

		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			curlF[1] = (-Rl.Fe*MVe - Rl.Fh*MVh)/2.0;
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			curlF[1] = (Rr.Fe*MVe + Rr.Fh*MVh)/2.0;
		}

	}

	
	// ====== Precompute terms ======

	std::complex<double> VIh, IVe;

	if (components_curl[2] || components_curl[3])
	{
		std::array<std::complex<double>, 10> terms_S1;
		ComputeTermsCurlKmi_SpatialS1(z, zp, rho, terms_S1);

		VIh = (Rl.Fh*terms_S1[1] + Rr.Fh*terms_S1[2] + Rl.Fh*Rr.Fh*(terms_S1[3] + terms_S1[4]) + Rm.Fh*(terms_S1[5] + Rl.Fh*terms_S1[6] + Rr.Fh*terms_S1[7] + Rl.Fh*Rr.Fh*(terms_S1[8] + terms_S1[9])))*MVh;

		IVe = (epsi/epsm)*(-Rl.Fe*terms_S1[1] - Rr.Fe*terms_S1[2] + Rl.Fe*Rr.Fe*(terms_S1[3] + terms_S1[4]) - Rm.Fe*(terms_S1[5] - Rl.Fe*terms_S1[6] - Rr.Fe*terms_S1[7] + Rl.Fe*Rr.Fe*(terms_S1[8] + terms_S1[9])))*MIe;
	}
	
	// For VV, the m < i case is treated as the reciprocal to the m > i case.
	
	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		reciprocal = true;
	}
	
	std::array<std::complex<double>, 10> terms_S0, terms_S2;
	std::complex<double> VVe_S2, VVe_S0, VVh_S2, VVh_S0;
	
	if (!reciprocal)
	{

		auto ComputeVVe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fe*terms[1] + Rr.Fe*terms[2] - Rl.Fe*Rr.Fe*(terms[3] - terms[4]) + this->Rm.Fe*(this->P*terms[5] - Rl.Fe*terms[6] + Rr.Fe*terms[7] - Rl.Fe*Rr.Fe*(terms[8] - terms[9])))*MVe;
			};

		auto ComputeVVh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fh*terms[1] + Rr.Fh*terms[2] - Rl.Fh*Rr.Fh*(terms[3] - terms[4]) + this->Rm.Fh*(this->P*terms[5] - Rl.Fh*terms[6] + Rr.Fh*terms[7] - Rl.Fh*Rr.Fh*(terms[8] - terms[9])))*MVh;
			};

		if (components_curl[0])
		{
			ComputeTermsCurlKmi_SpatialS2(z, zp, rho, terms_S2);
			VVe_S2 = ComputeVVe(terms_S2);
			VVh_S2 = ComputeVVh(terms_S2);
		}

		if (components_curl[1])
		{
			ComputeTermsCurlKmi_SpatialS0(z, zp, rho, terms_S0);
			VVe_S0 = ComputeVVe(terms_S0);
			VVh_S0 = ComputeVVh(terms_S0);
		}
		
	}
	else
	{

		auto ComputeIIe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((epsm/epsi)*(Rl.Fe*terms[1] + Rr.Fe*terms[2] + Rl.Fe*Rr.Fe*(terms[3] + terms[4]) - this->Rm.Fe*(terms[5] + Rl.Fe*terms[6] + Rr.Fe*terms[7] + Rl.Fe*Rr.Fe*(terms[8] + terms[9]))))*MVe;
			};

		auto ComputeIIh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((mui/mum)*(Rl.Fh*terms[1] + Rr.Fh*terms[2] + Rl.Fh*Rr.Fh*(terms[3] + terms[4]) - this->Rm.Fh*(terms[5] + Rl.Fh*terms[6] + Rr.Fh*terms[7] + Rl.Fh*Rr.Fh*(terms[8] + terms[9]))))*MVh;
			};

		if (components_curl[0])
		{
			ComputeTermsCurlKmi_SpatialS2(z, zp, rho, terms_S2);
			VVe_S2 = -ComputeIIe(terms_S2);
			VVh_S2 = -ComputeIIh(terms_S2);
		}

		if (components_curl[1])
		{
			ComputeTermsCurlKmi_SpatialS0(z, zp, rho, terms_S0);
			VVe_S0 = -ComputeIIe(terms_S0);
			VVh_S0 = -ComputeIIh(terms_S0);
		}
		
		// Reset
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		reciprocal = false;
		
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = VVe_S2 - VVh_S2;
	}

	if (components_curl[1])
	{
		curlK[1] = VVe_S0 + VVh_S0;
	}

	if (components_curl[2])
	{
		curlK[2] = VIh;
	}

	if (components_curl[3])
	{
		curlK[3] = IVe;
	}

	return;

}


/*! \brief Compute the GHJ-curl of the quasistatic spatial MGF when the source and observation layers are not the same.*/
void QuasistaticMGF::ComputeCurlGHJmi_Spatial(double z, double zp, double rho, std::array<std::complex<double>, 4> &curlK)
{	
	
	Fresnel Rr = GetFresnelCoefficient_Upward(i-1);
	Fresnel Rl = GetFresnelCoefficient_Downward(i);


	// ====== Compute singularity factors ======

	// Note: it is assumed that the homogeneous terms is always extracted
	if (extract_singularities)
	{

		std::fill(curlF.begin(), curlF.end(), 0.0);

		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			curlF[1] = (Rl.Fe*MVe + Rl.Fh*MVh)/2.0;
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			curlF[1] = (-Rr.Fe*MVe - Rr.Fh*MVh)/2.0;
		}

	}

	
	// ====== Precompute terms ======

	std::complex<double> VIh, IVe;

	if (components_curl[2] || components_curl[3])
	{
		std::array<std::complex<double>, 10> terms_S1;
		ComputeTermsCurlKmi_SpatialS1(z, zp, rho, terms_S1);

		VIh = (Rl.Fh*terms_S1[1] + Rr.Fh*terms_S1[2] + Rl.Fh*Rr.Fh*(terms_S1[3] + terms_S1[4]) + Rm.Fh*(terms_S1[5] + Rl.Fh*terms_S1[6] + Rr.Fh*terms_S1[7] + Rl.Fh*Rr.Fh*(terms_S1[8] + terms_S1[9])))*MVh;

		IVe = (epsi/epsm)*(-Rl.Fe*terms_S1[1] - Rr.Fe*terms_S1[2] + Rl.Fe*Rr.Fe*(terms_S1[3] + terms_S1[4]) - Rm.Fe*(terms_S1[5] - Rl.Fe*terms_S1[6] - Rr.Fe*terms_S1[7] + Rl.Fe*Rr.Fe*(terms_S1[8] + terms_S1[9])))*MIe;
	}
	
	// For VV, the m < i case is treated as the reciprocal to the m > i case.
	
	bool reciprocal = false;
	if (m < i)
	{
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		reciprocal = true;
	}
	
	std::array<std::complex<double>, 10> terms_S0, terms_S2;
	std::complex<double> IIe_S2, IIe_S0, IIh_S2, IIh_S0;
	
	if (!reciprocal)
	{

		auto ComputeIIe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((epsm/epsi)*(Rl.Fe*terms[1] + Rr.Fe*terms[2] + Rl.Fe*Rr.Fe*(terms[3] + terms[4]) - this->Rm.Fe*(terms[5] + Rl.Fe*terms[6] + Rr.Fe*terms[7] + Rl.Fe*Rr.Fe*(terms[8] + terms[9]))))*MVe;
			};

		auto ComputeIIh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (-1.0/2.0)*((mui/mum)*(Rl.Fh*terms[1] + Rr.Fh*terms[2] + Rl.Fh*Rr.Fh*(terms[3] + terms[4]) - this->Rm.Fh*(terms[5] + Rl.Fh*terms[6] + Rr.Fh*terms[7] + Rl.Fh*Rr.Fh*(terms[8] + terms[9]))))*MVh;
			};

		if (components_curl[0])
		{
			ComputeTermsCurlKmi_SpatialS2(z, zp, rho, terms_S2);
			IIe_S2 = ComputeIIe(terms_S2);
			IIh_S2 = ComputeIIh(terms_S2);
		}

		if (components_curl[1])
		{
			ComputeTermsCurlKmi_SpatialS0(z, zp, rho, terms_S0);
			IIe_S0 = ComputeIIe(terms_S0);
			IIh_S0 = ComputeIIh(terms_S0);
		}
		
	}
	else
	{

		auto ComputeVVe = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fe*terms[1] + Rr.Fe*terms[2] - Rl.Fe*Rr.Fe*(terms[3] - terms[4]) + this->Rm.Fe*(this->P*terms[5] - Rl.Fe*terms[6] + Rr.Fe*terms[7] - Rl.Fe*Rr.Fe*(terms[8] - terms[9])))*MVe;
			};

		auto ComputeVVh = [this, &Rr, &Rl] (std::array<std::complex<double>, 10> &terms)
			{
				return (1.0/2.0)*(-Rl.Fh*terms[1] + Rr.Fh*terms[2] - Rl.Fh*Rr.Fh*(terms[3] - terms[4]) + this->Rm.Fh*(this->P*terms[5] - Rl.Fh*terms[6] + Rr.Fh*terms[7] - Rl.Fh*Rr.Fh*(terms[8] - terms[9])))*MVh;
			};

		if (components_curl[0])
		{
			ComputeTermsCurlKmi_SpatialS2(z, zp, rho, terms_S2);
			IIe_S2 = -ComputeVVe(terms_S2);
			IIh_S2 = -ComputeVVh(terms_S2);
		}

		if (components_curl[1])
		{
			ComputeTermsCurlKmi_SpatialS0(z, zp, rho, terms_S0);
			IIe_S0 = -ComputeVVe(terms_S0);
			IIh_S0 = -ComputeVVh(terms_S0);
		}
		
		// Reset
		SwapSourceAndObserver(z, zp);
		Rr = GetFresnelCoefficient_Upward(i-1);
		Rl = GetFresnelCoefficient_Downward(i);
		reciprocal = false;
		
	}

	
	// ====== Assemble required component(s) ======

	if (components_curl[0])
	{
		curlK[0] = IIh_S2 - IIe_S2;
	}

	if (components_curl[1])
	{
		curlK[1] = IIh_S0 + IIe_S0;
	}

	if (components_curl[2])
	{
		curlK[2] = IVe;
	}

	if (components_curl[3])
	{
		curlK[3] = VIh;
	}

	return;

}


// ==================================================================================
// Computational helpers
// ==================================================================================

/*! \brief Compute the z- and zp-dependent same-layer spectral terms for the 0-order Sommerfeld integral.*/
void QuasistaticMGF::ComputeExpTermsKii_kz(double z, double zp, std::complex<double> kzi, std::array<std::complex<double>, 5> &exp_terms)
{

	std::complex<double> jkzi = J*kzi;

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	for (int ii = 0; ii < exp_terms.size(); ii++)
		exp_terms[ii] = std::exp(-jkzi*gamma[ii]);
	
	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spectral terms for the 1st-order Sommerfeld integral.*/
void QuasistaticMGF::ComputeExpTermsKii_krho(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 5> &exp_terms)
{

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	for (int ii = 0; ii < exp_terms.size(); ii++)
		exp_terms[ii] = std::exp(-krho*gamma[ii]);

	return;

}


/*! \brief Compute the z- and zp-dependent different-layer spectral terms, for the krho-form of all terms.*/
void QuasistaticMGF::ComputeExpTermsKmi_krho(double z, double zp, std::complex<double> krho, std::array<std::complex<double>, 10> &exp_terms)
{

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);
	
	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	for (int ii = 0; ii < alpha.size(); ii++)
		exp_terms[ii] = std::exp(-krho*alpha[ii]);

	return;

}


/*! \brief Compute the z- and zp-dependent different-layer spectral terms, for the kz-form for S0 terms.*/
void QuasistaticMGF::ComputeExpTermsKmi_kz(double z, double zp, std::complex<double> kzi, std::array<std::complex<double>, 10> &exp_terms)
{

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);
		
	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	for (int ii = 0; ii < alpha.size(); ii++)
		exp_terms[ii] = std::exp(-J*kzi*alpha[ii]);

	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 0-order Sommerfeld integral.*/
void QuasistaticMGF::ComputeExpTermsKii_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 5> &exp_terms)
{

	std::complex<double> jki = J*ki;

	std::array<double, 5> ksi;
	ComputeDistanceTermsKii(z, zp, rho, ksi);
	
	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);

	if (!extract_singularities)
	{
		for (int ii = 0; ii < exp_terms.size(); ii++)
			exp_terms[ii] = std::exp(-jki*ksi[ii])/ksi[ii];
	}
	else
	{
		exp_terms[0] = ComputeNonsingularHGF(ksi[0], ki);

		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			exp_terms[1] = exp_terms[0];
			exp_terms[2] = std::exp(-jki*ksi[2])/ksi[2];
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			exp_terms[1] = std::exp(-jki*ksi[1])/ksi[1];
			exp_terms[2] = exp_terms[0];
		}
		else
		{
			exp_terms[1] = std::exp(-jki*ksi[1])/ksi[1];
			exp_terms[2] = std::exp(-jki*ksi[2])/ksi[2];
		}

		exp_terms[3] = std::exp(-jki*ksi[3])/ksi[3];
		exp_terms[4] = std::exp(-jki*ksi[4])/ksi[4];
	}
	
	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 1st-order Sommerfeld integral.*/
void QuasistaticMGF::ComputeExpTermsKii_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 5> &exp_terms)
{

	std::array<double, 5> ksi;
	ComputeDistanceTermsKii(z, zp, rho, ksi);

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	
	// If rho ~ 0, J_1(rho*krho) ~ 0, so off-diagonal terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	if (!extract_singularities)
	{
		for (int ii = 1; ii < exp_terms.size(); ii++)
			exp_terms[ii] = (1.0/rho)*(1.0 - gamma[ii]/ksi[ii]);
	}
	else
	{
		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			exp_terms[1] = 0.0;
			exp_terms[2] = (1.0/rho)*(1.0 - gamma[2]/ksi[2]);
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			exp_terms[1] = (1.0/rho)*(1.0 - gamma[1]/ksi[1]);
			exp_terms[2] = 0.0;
		}
		else
		{
			exp_terms[1] = (1.0/rho)*(1.0 - gamma[1]/ksi[1]);
			exp_terms[2] = (1.0/rho)*(1.0 - gamma[2]/ksi[2]);
		}
		
		exp_terms[3] = (1.0/rho)*(1.0 - gamma[3]/ksi[3]);
		exp_terms[4] = (1.0/rho)*(1.0 - gamma[4]/ksi[4]);
	}

	return;

}


// /*! \brief Compute the z- and zp-dependent different-layer spatial terms for the 0-order Sommerfeld integral. If using krho-form of the Sommerfeld integrand.*/
// void QuasistaticMGF::ComputeExpTermsKmi_SpatialS0(double z, double zp, double rho, std::vector<std::complex<double>> &exp_terms)
// {

// 	std::vector<double> ksi (10);
// 	ComputeDistanceTermsKmi(z, zp, rho, ksi);
	
// 	exp_terms.resize(10);

// 	if (!extract_singularities || std::abs(m - i) > 1)
// 	{
// 		for (int ii = 0; ii < exp_terms.size(); ii++)
// 			exp_terms[ii] = 1.0/ksi[ii];
// 	}
// 	else if (std::abs(m - i) == 1)
// 	{
// 		exp_terms[0] = 0.0;
		
// 		if (!AxiallyCoincident(z, zp))
// 		{
// 			for (int ii = 1; ii < exp_terms.size(); ii++)
// 				exp_terms[ii] = 1.0/ksi[ii];
// 		}
// 		else
// 		{
// 			if (m == i + 1)
// 			{
// 				exp_terms[1] = 0.0;
// 				exp_terms[2] = 1.0/ksi[2];
// 			}
// 			else if (m == i - 1)
// 			{
// 				exp_terms[1] = 1.0/ksi[1];
// 				exp_terms[2] = 0.0;
// 			}

// 			for (int ii = 3; ii < exp_terms.size(); ii++)
// 				exp_terms[ii] = 1.0/ksi[ii];
// 		}
// 	}
	
// 	return;

// }


/*! \brief Compute the z- and zp-dependent different-layer spatial terms for the 0-order Sommerfeld integral. If using kz-form of the Sommerfeld integrand.*/
void QuasistaticMGF::ComputeExpTermsKmi_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 10> &exp_terms)
{

	std::complex<double> jki = J*ki;
	
	std::array<double, 10> ksi;
	ComputeDistanceTermsKmi(z, zp, rho, ksi);
	
	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);

	if (!extract_singularities)
	{
		for (int ii = 0; ii < exp_terms.size(); ii++)
			exp_terms[ii] = std::exp(-jki*ksi[ii])/ksi[ii];
	}
	else
	{
		exp_terms[0] = ComputeNonsingularHGF(ksi[0], ki);
		
		if (!AxiallyCoincident(z, zp))
		{
			for (int ii = 1; ii < exp_terms.size(); ii++)
				exp_terms[ii] = std::exp(-jki*ksi[ii])/ksi[ii];
		}
		else
		{
			if (m == i + 1)
			{
				exp_terms[1] = exp_terms[0];
				exp_terms[2] = std::exp(-jki*ksi[2])/ksi[2];;
			}
			else if (m == i - 1)
			{
				exp_terms[1] = std::exp(-jki*ksi[1])/ksi[1];;
				exp_terms[2] = exp_terms[0];
			}

			for (int ii = 3; ii < exp_terms.size(); ii++)
				exp_terms[ii] = std::exp(-jki*ksi[ii])/ksi[ii];
		}
	}
	
	return;

}


/*! \brief Compute the z- and zp-dependent different-layer spatial terms for the 1st-order Sommerfeld integral.*/
void QuasistaticMGF::ComputeExpTermsKmi_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 10> &exp_terms)
{
	
	std::array<double, 10> ksi;
	ComputeDistanceTermsKmi(z, zp, rho, ksi);

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);
	
	std::fill(exp_terms.begin(), exp_terms.end(), 0.0);
	
	// If rho ~ 0, J_1(rho*krho) ~ 0, so off-diagonal terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	if (!extract_singularities || !AxiallyCoincident(z, zp))
	{
		for (int ii = 0; ii < exp_terms.size(); ii++)
			exp_terms[ii] = (1.0/rho)*(1.0 - alpha[ii]/ksi[ii]);
	}
	else if (AxiallyCoincident(z, zp))
	{
		exp_terms[0] = 0.0;
		
		if (m == i + 1)
		{
			exp_terms[1] = 0.0;
			exp_terms[2] = (1.0/rho)*(1.0 - alpha[2]/ksi[2]);
		}
		else if (m == i - 1)
		{
			exp_terms[1] = (1.0/rho)*(1.0 - alpha[1]/ksi[1]);
			exp_terms[2] = 0.0;
		}

		for (int ii = 3; ii < exp_terms.size(); ii++)
			exp_terms[ii] = (1.0/rho)*(1.0 - alpha[ii]/ksi[ii]);
	}

	return;

}


/*! \brief Helper function to compute same-layer distance terms.*/
void QuasistaticMGF::ComputeDistanceTermsKii(double z, double zp, double rho, std::array<double, 5> &ksi)
{
	
	auto ksi_fn = [&rho](double _gamma) -> double
		{
			return std::sqrt(std::pow(rho, 2) + std::pow(_gamma, 2));
		};

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(ksi.begin(), ksi.end(), 0.0);
	for (int ii = 0; ii < ksi.size(); ii++)
		ksi[ii] = ksi_fn(gamma[ii]);

	return;

}


/*! \brief Helper function to compute different-layer distance terms.*/
void QuasistaticMGF::ComputeDistanceTermsKmi(double z, double zp, double rho, std::array<double, 10> &ksi)
{

	auto ksi_fn = [&rho](double _alpha) -> double
		{
			return std::sqrt(std::pow(rho, 2) + std::pow(_alpha, 2));
		};

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);

	std::fill(ksi.begin(), ksi.end(), 0.0);
	for (int ii = 0; ii < ksi.size(); ii++)
		ksi[ii] = ksi_fn(alpha[ii]);

	return;

}


/*! \brief Helper function to compute same-layer z-zp terms.*/
void QuasistaticMGF::ComputeGammaTerms(double z, double zp, std::array<double, 5> &gamma)
{

	// double zmin = lm->layers[i].zmin;
	// double zmax = lm->layers[i].zmax;
	// double h = lm->layers[i].h;
	double zmin = lm->GetZmin(i);
	double zmax = lm->GetZmax(i);
	double h = lm->GetHeight(i);
	
	gamma[0] = std::abs(z - zp);
	gamma[1] = z + zp - 2.0*zmin;
	gamma[2] = 2.0*zmax - (z + zp);
	gamma[3] = 2.0*h - z + zp;
	gamma[4] = 2.0*h + z - zp;

	return;

}


/*! \brief Helper function to compute different-layer z terms.*/
void QuasistaticMGF::ComputeTauTerms(double z, std::array<double, 2> &tau)
{

	// double zmin = lm->layers[m].zmin;
	// double zmax = lm->layers[m].zmax;
	double zmin = lm->GetZmin(m);
	double zmax = lm->GetZmax(m);

	std::fill(tau.begin(), tau.end(), 0.0);

	if (m < i)
	{
		tau[0] = z - zmin;
		tau[1] = -z - zmin + 2.0*zmax;
	}
	else if (m > i)
	{
		tau[0] = zmax - z;
		tau[1] = z - 2.0*zmin + zmax;
	}

	return;

}


/*! \brief Helper function to compute different-layer z-zp terms.*/
void QuasistaticMGF::ComputeAlphaTerms(double z, double zp, std::array<double, 10> &alpha)
{

	// double zmin = lm->layers[i].zmin;
	// double zmax = lm->layers[i].zmax;
	double zmin = lm->GetZmin(i);
	double zmax = lm->GetZmax(i);

	double _z;
	if (m < i)
		_z = zmax;
	else if (m > i)
		_z = zmin;

	std::array<double, 5> gamma;
	ComputeGammaTerms(_z, zp, gamma);

	std::array<double, 2> tau;
	ComputeTauTerms(z, tau);
	
	std::fill(alpha.begin(), alpha.end(), 0.0);
	for (int idx_tau = 0; idx_tau < tau.size(); idx_tau++)
		for (int idx_gamma = 0; idx_gamma < gamma.size(); idx_gamma++)
			alpha[idx_gamma + 5*idx_tau] = gamma[idx_gamma] + tau[idx_tau] + nu;
	
	return;

}


/*! \brief Compute the upwards reflection coefficient multipliers, between the observation layer m and the source layer i of the equivalent transmission line network.*/
void QuasistaticMGF::ComputeM_Upward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh, double &_nu)
{

	_MVe = 1.0;
	_MVh = 1.0;
	_MIe = 1.0;
	_MIh = 1.0;

	_nu = 0.0;
	
	for (int ii = m; ii <= i-2; ii++)
	{
		Fresnel Rr = GetFresnelCoefficient_Upward(ii);
		
		_MVe *= (1.0 + Rr.Fe);
		_MVh *= (1.0 + Rr.Fh);

		_MIe *= (1.0 - Rr.Fe);
		_MIh *= (1.0 - Rr.Fh);

		// _nu += lm->layers[ii+1].h;
		_nu += lm->GetHeight(ii+1);
	}

	return;

}


/*! \brief Compute the downwards reflection coefficient multipliers, between the observation layer m and the source layer i of the equivalent transmission line network.*/
void QuasistaticMGF::ComputeM_Downward(std::complex<double> &_MVe, std::complex<double> &_MVh, std::complex<double> &_MIe, std::complex<double> &_MIh, double &_nu)
{

	_MVe = 1.0;
	_MVh = 1.0;
	_MIe = 1.0;
	_MIh = 1.0;

	_nu = 0.0;
	
	for (int ii = i+1; ii <= m-1; ii++)
	{
		Fresnel Rl = GetFresnelCoefficient_Downward(ii);
			     
		_MVe *= (1.0 + Rl.Fe);
		_MVh *= (1.0 + Rl.Fh);

		_MIe *= (1.0 - Rl.Fe);
		_MIh *= (1.0 - Rl.Fh);

		// _nu += lm->layers[ii].h;
		_nu += lm->GetHeight(ii);
	}

	return;

}


/*! \brief Function to compute the nonsingular part of the homogeneous Green's function, exp(-jkr)/r - 1/r, using a Taylor expansion to remove the singularity.*/
std::complex<double> QuasistaticMGF::ComputeNonsingularHGF(double r, std::complex<double> k)
{

	if (std::abs(k) < 1.0e-15)
		return 0.0;
	else if (std::abs(k*r) > 0.1)
		return (std::exp(-J*k*r) - 1.0)/r;

	std::complex<double> result;
	double threshold = 1.0e-8;
	// int max_it = 40;

	std::complex<double> prev = -J*k;
	result = prev;

	for (int ii = 2; true; ii++)
	{
		prev = prev*(-J*k*r)/(double)ii;
		result += prev;
			
		// if (std::abs(prev)/std::abs(result) < threshold || ii > max_it)
		if (std::abs(prev)/std::abs(result) < threshold)
			break;
	}
	
	return result;
	
}


/*! \brief Function to test whether the source and observation z-points are the same, or close enough to be considered the same.*/
bool QuasistaticMGF::AxiallyCoincident(double z, double zp, double tol)
{

	// Normalize the comparison by the total height of the stackup
	tol *= std::abs(lm->layers[0].zmax - lm->layers.back().zmin);
	if (std::abs(z - zp) > tol)
		return false;
	else
		return true;

}


/*! \brief Function to test whether the source and observation x- and y-points are the same, or close enough to be considered the same.*/
bool QuasistaticMGF::LaterallyCoincident(double rho, double tol)
{

	// Normalize the comparison by the total height of the stackup
	tol *= std::abs(lm->layers[0].zmax - lm->layers.back().zmin);
	if (std::abs(rho) > tol)
		return false;
	else
		return true;

}


// ==================================================================================
// Computational helpers - curl
// ==================================================================================

/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 0-order Sommerfeld integral, for the curl components.*/
void QuasistaticMGF::ComputeTermsCurlKii_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms)
{

	std::complex<double> jki = J*ki;
	
	std::array<double, 5> ksi;
	ComputeDistanceTermsKii(z, zp, rho, ksi);

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(terms.begin(), terms.end(), 0.0);
	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = gamma[ii]*(1.0 + jki*ksi[ii])*std::exp(-jki*ksi[ii])/std::pow(ksi[ii], 3);

	if (extract_singularities && AxiallyCoincident(z, zp))
	{
		// Source and observer are at the bottom interface of the layer
		// if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmin))
		if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmin(i)))
		{
			terms[1] = 0.0;
		}
		// Source and observer are at the top interface of the layer
		// else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->layers[i].zmax))
		else if (AxiallyCoincident(z, zp) && AxiallyCoincident(z, lm->GetZmax(i)))
		{
			terms[2] = 0.0;
		}
	}
	
	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 1st-order Sommerfeld integral, for the curl components.*/
void QuasistaticMGF::ComputeTermsCurlKii_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms)
{

	std::complex<double> jki = J*ki;
	
	std::array<double, 5> ksi;
	ComputeDistanceTermsKii(z, zp, rho, ksi);

	std::fill(terms.begin(), terms.end(), 0.0);
	
	// If rho ~ 0, J_1(rho*krho) ~ 0, so spatial terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = rho*(1.0 + jki*ksi[ii])*std::exp(-jki*ksi[ii])/std::pow(ksi[ii], 3);

	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 1st-order Sommerfeld integral, for the curl components.*/
void QuasistaticMGF::ComputeTermsCurlKii_SpatialS2(double z, double zp, double rho, std::array<std::complex<double>, 5> &terms)
{

	std::array<double, 5> ksi;
	ComputeDistanceTermsKii(z, zp, rho, ksi);

	std::array<double, 5> gamma;
	ComputeGammaTerms(z, zp, gamma);

	std::fill(terms.begin(), terms.end(), 0.0);
	
	// If rho ~ 0, J_2(rho*krho) ~ 0, so spatial terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = std::pow((-gamma[ii] + ksi[ii]), 2)*(2.0*ksi[ii] + gamma[ii])/std::pow(rho, 2)/std::pow(ksi[ii], 3);

	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 0-order Sommerfeld integral, for different-layer curl components.*/
void QuasistaticMGF::ComputeTermsCurlKmi_SpatialS0(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms)
{

	std::complex<double> jki = J*ki;

	std::array<double, 10> ksi;
	ComputeDistanceTermsKmi(z, zp, rho, ksi);

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);

	std::fill(terms.begin(), terms.end(), 0.0);
	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = alpha[ii]*(1.0 + jki*ksi[ii])*std::exp(-jki*ksi[ii])/std::pow(ksi[ii], 3);

	if (extract_singularities && AxiallyCoincident(z, zp))
	{
		if (m == i + 1)
		{
			terms[1] = 0.0;
		}
		else if (m == i - 1)
		{
			terms[2] = 0.0;
		}
	}
	
	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 1st-order Sommerfeld integral, for different-layer curl components.*/
void QuasistaticMGF::ComputeTermsCurlKmi_SpatialS1(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms)
{

	std::complex<double> jki = J*ki;

	std::array<double, 10> ksi;
	ComputeDistanceTermsKmi(z, zp, rho, ksi);
	
	std::fill(terms.begin(), terms.end(), 0.0);
	
	// If rho ~ 0, J_1(rho*krho) ~ 0, so spatial terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = rho*(1.0 + jki*ksi[ii])*std::exp(-jki*ksi[ii])/std::pow(ksi[ii], 3);

	return;

}


/*! \brief Compute the z- and zp-dependent same-layer spatial terms for the 1st-order Sommerfeld integral, for different-layer curl components.*/
void QuasistaticMGF::ComputeTermsCurlKmi_SpatialS2(double z, double zp, double rho, std::array<std::complex<double>, 10> &terms)
{

	std::array<double, 10> ksi;
	ComputeDistanceTermsKmi(z, zp, rho, ksi);

	std::array<double, 10> alpha;
	ComputeAlphaTerms(z, zp, alpha);
	
	std::fill(terms.begin(), terms.end(), 0.0);
	
	// If rho ~ 0, J_2(rho*krho) ~ 0, so spatial terms are ~0.
	if (LaterallyCoincident(rho))
		return;

	for (int ii = 1; ii < terms.size(); ii++)
		terms[ii] = std::pow((-alpha[ii] + ksi[ii]), 2)*(2.0*ksi[ii] + alpha[ii])/std::pow(rho, 2)/std::pow(ksi[ii], 3);

	return;

}


/*! \brief Swap the source and observation point and redo all associated precomputations. Useful when using reciprocity to compute components.*/
void QuasistaticMGF::SwapSourceAndObserver(double &z, double &zp)
{

	int i_old = i, m_old = m;
	double z_old = z, zp_old = zp;
	SetLayers(m_old, i_old);
	zp = z_old;
	z = zp_old;

	if (z > zp)
		P = 1.0;
	else if (z < zp)
		P = -1.0;
	else
		P = 0.0;
	
	return;
	
}



