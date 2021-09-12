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

/********************************** MGF.cpp *******************************

 * Routines for computing the multilayer Green's function (MGF) based on
 * Michalski & Zheng's Formulation-C. Currently supports direct numerical
 * integration, multilevel DCIM, interpolation table usage and re-usage, 
 * and quasistatic analysis and extraction.
 *
 * Author: Shashwat Sharma
 * Created on: Apr 02, 2020

 *************************************************************************/


#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>

#include "MGF.hpp"
#include "layers.hpp"
#include "spectral_MGF.hpp"
#include "quasistatic_MGF.hpp"
#include "sommerfeld_integrals.hpp"
#include "DCIM.hpp"
#include "constants.hpp"

using namespace strata;


// ==================================================================================
// Interface
// ==================================================================================

/*! \brief Initialize the MGF engine for the given frequency and settings.*/
void MGF::Initialize(double _f, LayerManager &_lm, MGF_settings &_s)
{
	
	s = _s;

	if (s.verbose)
		std::cout << "\n> Initializing MGF..." << std::endl;
	
	f = _f;
	lm = _lm;
	omega = 2.0*M_PI*f;

	lm.ProcessLayers(f);


	// ====== Initialize quasistatic MGF ======
	
	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
	{
		qmgf.Initialize(&lm, f, s.components, s.extract_singularities);

		if (s.curl_type == MGF_GEM)
			qmgf.SetCurlToGEM();
		else
			qmgf.SetCurlToGHJ();
	}


	// ====== Initialize spectral MGF ======

	if (s.method == MGF_INTEGRATE || s.method == MGF_INTERPOLATE)
	{
		smgf.Initialize(&lm, f, s.components, s.extract_quasistatic, s.extract_homogeneous);
		
		if (s.switching_point < 0.0)
			s.switching_point = 1.2*std::abs(lm.k_max);

		if (s.curl_type == MGF_GEM)
			smgf.SetCurlToGEM();
		else
			smgf.SetCurlToGHJ();
	}
	
	
	// ====== Generate DCIM images ======

	if (s.method == MGF_DCIM ||
		(s.method == MGF_INTERPOLATE && s.sampling_method == MGF_DCIM))
	{
		if (s.verbose)
			std::cout << "--- Generating DCIM images..." << std::flush;
		
		DCIM_settings s_dcim;
		s_dcim.components = s.components;
		s_dcim.method = s.DCIM_method;
		s_dcim.extract_quasistatic = s.extract_quasistatic;
		s_dcim.extract_homogeneous = s.extract_homogeneous;

		s_dcim.components_curl = s.components_curl;
		s_dcim.compute_curl = s.compute_curl;
		if (s.curl_type == MGF_GEM)
			s_dcim.curl_GEM = true;
		else if (s.curl_type == MGF_GHJ)
			s_dcim.curl_GEM = false;
		
		s_dcim.tol_svd = s.tol_svd;
		s_dcim.tol_eig = s.tol_eig;
		s_dcim.max_num_images = s.max_num_images;

		dcim.GenerateImages(&lm, f, s_dcim);

		if (s.verbose)
			std::cout << "done." << std::endl;
	}

	
	// ====== Build or import interpolation table ======

	if (s.method == MGF_INTERPOLATE)
	{

		// ====== MGF ======
		
		int load_info = -1, export_info = -1;
		
		if (s.load_table && !s.export_table)
		{
			if (s.verbose)
				std::cout << "--- Loading interpolation table..." << std::flush;

			load_info = LoadTable(MGF_table, s.filename);

			if (s.verbose && load_info == 0)
				std::cout << "done." << std::endl;

		}
		
		if (!s.load_table || load_info != 0 || s.export_table)
		{
			if (s.verbose)
				std::cout << "--- Generating interpolation table..." << std::flush;

			TabulateMGF(MGF_table, false);

			if (s.verbose)
				std::cout << "done." << std::endl;

			if (s.export_table)
			{
				if (s.verbose)
					std::cout << "--- Exporting interpolation table..." << std::flush;

				export_info = ExportTable(MGF_table, s.filename);

				if (s.verbose)
					std::cout << "done." << std::endl;
			}
		}


		// ====== MGF curl ======

		if (s.compute_curl)
		{
			
			int load_info_curl = -1, export_info_curl = -1;
		
			if (s.load_table && !s.export_table)
			{
				if (s.verbose)
					std::cout << "--- Loading curl interpolation table..." << std::flush;

				load_info_curl = LoadTable(CurlMGF_table, s.filename_curl);

				if (s.verbose && load_info_curl == 0)
					std::cout << "done." << std::endl;

			}
		
			if (!s.load_table || load_info_curl != 0 || s.export_table)
			{
				if (s.verbose)
					std::cout << "--- Generating curl interpolation table..." << std::flush;

				TabulateMGF(CurlMGF_table, true);

				if (s.verbose)
					std::cout << "done." << std::endl;

				if (s.export_table)
				{
					if (s.verbose)
						std::cout << "--- Exporting curl interpolation table..." << std::flush;

					export_info_curl = ExportTable(CurlMGF_table, s.filename_curl);

					if (s.verbose)
						std::cout << "done." << std::endl;
				}
			}

		}
		
	}


	// ====== Update switchboard ======
	
	initialized = true;
	layers_set = false;
	singularity_factors_computed = false;

	if (s.verbose)
		std::cout << "> done.\n" << std::endl;

	return;
	
}


/*! \brief Set source and observation layers and execute related precomputations.*/
void MGF::SetLayers(int _i, int _m)
{

	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF::SetLayers(): Must call MGF::Initialize() before setting layer indices.");
	}

	i = _i;
	m = _m;

	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
		qmgf.SetLayers(i, m);

	if (s.method == MGF_INTEGRATE)
		smgf.SetLayers(i, m);

	layers_set = true;
	singularity_factors_computed = false;
	
	return;
	
}


/*! \brief Update singularity extraction settings.*/
void MGF::SetSingularityExtraction(bool extract_singularities)
{
	
	s.extract_singularities = extract_singularities;
	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
		qmgf.extract_singularities = extract_singularities;
	
	return;
	
}


/*! \brief Update components to compute.*/
void MGF::SetComponents(std::vector<bool> components)
{
	
	s.components = components;
	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
		qmgf.components = components;
	if (s.method == MGF_INTEGRATE || s.method == MGF_INTERPOLATE)
		smgf.components = components;
	if (s.method == MGF_DCIM ||
		(s.method == MGF_INTERPOLATE && s.sampling_method == MGF_DCIM))
		dcim.s.components = components;
	
	return;
	
}


/*! \brief Update curl components to compute.*/
void MGF::SetCurlComponents(std::vector<bool> components_curl)
{
	
	s.components_curl = components_curl;
	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
		qmgf.components_curl = components_curl;
	if (s.method == MGF_INTEGRATE || s.method == MGF_INTERPOLATE)
		smgf.components_curl = components_curl;
	if (s.method == MGF_DCIM ||
		(s.method == MGF_INTERPOLATE && s.sampling_method == MGF_DCIM))
		dcim.s.components_curl = components_curl;
	
	return;
	
}


/*! \brief Compute the MGF. Note that x_diff = x_obs - x_src, and y_diff = y_obs - y_src. MGF components are stored in row-major order.*/
void MGF::ComputeMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G, std::complex<double> &G_phi)
{

	// ====== Sanity checks ======
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF::ComputeMGF(): Must call MGF::Initialize() before computing the MGF.");
	}

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] MGF::ComputeMGF(): Must call MGF::SetLayers() before computing the MGF.");
	}

	// if ((z > lm.layers[0].zmax || z < lm.layers.back().zmin) ||
	// 	(zp > lm.layers[0].zmax || zp < lm.layers.back().zmin))
	// {
	// 	throw std::domain_error("[ERROR] MGF::ComputeMGF(): Source or observation z-point is out of bounds.");
	// }


	// ====== Computation ======
	
	double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

	std::complex<double> zeta = std::atan2(y_diff, x_diff);
	cos_term = std::cos(zeta);
	sin_term = std::sin(zeta);

	std::array<std::complex<double>, 5> _G;
	std::fill(_G.begin(), _G.end(), 0.0);

	if (s.method == MGF_INTEGRATE)
	{
		ComputeMGF_Integration(rho, z, zp, _G);
	}
	else if (s.method == MGF_DCIM)
	{
		ComputeMGF_DCIM(rho, z, zp, _G);
	}
	else if (s.method == MGF_INTERPOLATE)
	{
		ComputeMGF_Interpolation(rho, z, zp, _G, MGF_table, s.components);
	}

	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
	{
		qmgf.ComputeQMGF_Spatial(z, zp, rho);

		for (int jj = 0; jj < _G.size(); jj++)
			_G[jj] += qmgf.GetResult_Spatial(jj);

		if (s.extract_singularities || std::abs(lm.k_max) < 1.0e-15)
			ComputeSingularityFactors();
	}
	else if (s.extract_homogeneous)
		ComputeHomogeneousFactors();

	G[0] = _G[0];
	G[1] = 0.0;
	G[2] = _G[1]*cos_term;
	G[3] = 0.0;
	G[4] = _G[0];
	G[5] = _G[1]*sin_term;
	G[6] = _G[2]*cos_term;
	G[7] = _G[2]*sin_term;
	G[8] = _G[3];
	
	G_phi = _G[4];
	
	return;

}


/*! \brief Compute the curl components of the MGF. Note that x_diff = x_obs - x_src, and y_diff = y_obs - y_src. MGF components are stored in row-major order.*/
void MGF::ComputeCurlMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G)
{

	// ====== Sanity checks ======
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF::ComputeCurlMGF(): Must call MGF::Initialize() before computing the MGF.");
	}

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] MGF::ComputeCurlMGF(): Must call MGF::SetLayers() before computing the MGF.");
	}

	// if ((z > lm.layers[0].zmax || z < lm.layers.back().zmin) ||
	// 	(zp > lm.layers[0].zmax || zp < lm.layers.back().zmin))
	// {
	// 	throw std::domain_error("[ERROR] MGF::ComputeCurlMGF(): Source or observation z-point is out of bounds.");
	// }


	// ====== Computation ======
	
	double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

	std::complex<double> zeta = std::atan2(y_diff, x_diff);
	cos_term = std::cos(zeta);
	sin_term = std::sin(zeta);
	cos2_term = std::cos(2.0*zeta);
	sin2_term = std::sin(2.0*zeta);

	std::array<std::complex<double>, 4> _G;
	std::fill(_G.begin(), _G.end(), 0.0);

	if (s.method == MGF_INTEGRATE)
	{
		ComputeCurlMGF_Integration(rho, z, zp, _G);
	}
	else if (s.method == MGF_DCIM)
	{
		ComputeCurlMGF_DCIM(rho, z, zp, _G);
	}
	else if (s.method == MGF_INTERPOLATE)
	{
		ComputeMGF_Interpolation(rho, z, zp, _G, CurlMGF_table, s.components_curl);
	}

	if (s.extract_quasistatic || s.method == MGF_QUASISTATIC)
	{
		qmgf.ComputeCurlQMGF_Spatial(z, zp, rho);

		for (int jj = 0; jj < _G.size(); jj++)
			_G[jj] += qmgf.GetResultCurl_Spatial(jj);
	}

	G[0] = -(_G[0]*sin2_term/2.0);
	G[1] = (_G[0]*cos2_term - _G[1])/2.0;
	G[2] = (_G[2]*sin_term);
	G[3] = (_G[0]*cos2_term + _G[1])/2.0;
	G[4] = (_G[0]*sin2_term)/2.0;
	G[5] = -(_G[2]*cos_term);
	G[6] = -(_G[3]*sin_term);
	G[7] = (_G[3]*cos_term);
	G[8] = 0.0;

	// for (std::size_t ii = 0; ii < G.size(); ii++)
	// 	G[ii] *= -1.0;
	
	return;

}


/*! \brief Compute the quasistatic MGF. Note that x_diff = x_obs - x_src, and y_diff = y_obs - y_src. MGF components are stored in row-major order.*/
void MGF::ComputeQMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G, std::complex<double> &G_phi)
{

	// ====== Sanity checks ======
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF::ComputeQMGF(): Must call MGF::Initialize() before computing the quasistatic MGF.");
	}

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] MGF::ComputeQMGF(): Must call MGF::SetLayers() before computing the quasistatic MGF.");
	}

	// if ((z > lm.layers[0].zmax || z < lm.layers.back().zmin) ||
	// 	(zp > lm.layers[0].zmax || zp < lm.layers.back().zmin))
	// {
	// 	throw std::domain_error("[ERROR] MGF::ComputeQMGF(): Source or observation z-point is out of bounds.");
	// }


	// ====== Computation ======
	
	double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

	cos_term = std::cos(std::atan2(y_diff, x_diff));
	sin_term = std::sin(std::atan2(y_diff, x_diff));

	std::array<std::complex<double>, 5> _G;
	std::fill(_G.begin(), _G.end(), 0.0);

	qmgf.ComputeQMGF_Spatial(z, zp, rho);

	for (int jj = 0; jj < _G.size(); jj++)
		_G[jj] += qmgf.GetResult_Spatial(jj);

	if (s.extract_singularities || std::abs(lm.k_max) < 1.0e-15)
		ComputeSingularityFactors();

	// _G[1] = 0.0; _G[2] = 0.0;
	
	G[0] = _G[0];
	G[1] = 0.0;
	G[2] = _G[1]*cos_term;
	G[3] = 0.0;
	G[4] = _G[0];
	G[5] = _G[1]*sin_term;
	G[6] = _G[2]*cos_term;
	G[7] = _G[2]*sin_term;
	G[8] = _G[3];
	
	G_phi = _G[4];
	
	return;

}


/*! \brief Compute the curl components of the quasistatic MGF. Note that x_diff = x_obs - x_src, and y_diff = y_obs - y_src. Quasistatic MGF components are stored in row-major order.*/
void MGF::ComputeCurlQMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G)
{

	// ====== Sanity checks ======
	
	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF::ComputeCurlMGF(): Must call MGF::Initialize() before computing the MGF.");
	}

	if (!layers_set)
	{
		throw std::logic_error("[ERROR] MGF::ComputeCurlMGF(): Must call MGF::SetLayers() before computing the MGF.");
	}

	// if ((z > lm.layers[0].zmax || z < lm.layers.back().zmin) ||
	// 	(zp > lm.layers[0].zmax || zp < lm.layers.back().zmin))
	// {
	// 	throw std::domain_error("[ERROR] MGF::ComputeCurlMGF(): Source or observation z-point is out of bounds.");
	// }


	// ====== Computation ======
	
	double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

	std::complex<double> zeta = std::atan2(y_diff, x_diff);
	cos_term = std::cos(zeta);
	sin_term = std::sin(zeta);
	cos2_term = std::cos(2.0*zeta);
	sin2_term = std::sin(2.0*zeta);

	std::array<std::complex<double>, 4> _G;
	std::fill(_G.begin(), _G.end(), 0.0);

	qmgf.ComputeCurlQMGF_Spatial(z, zp, rho);

	for (int jj = 0; jj < _G.size(); jj++)
		_G[jj] += qmgf.GetResultCurl_Spatial(jj);

	G[0] = -(_G[0]*sin2_term/2.0);
	G[1] = (_G[0]*cos2_term - _G[1])/2.0;
	G[2] = (_G[2]*sin_term);
	G[3] = (_G[0]*cos2_term + _G[1])/2.0;
	G[4] = (_G[0]*sin2_term)/2.0;
	G[5] = -(_G[2]*cos_term);
	G[6] = -(_G[3]*sin_term);
	G[7] = (_G[3]*cos_term);
	G[8] = 0.0;

	return;

}


/*! \brief Accessor function to return the singularity factor for a particular MGF component. Assumes ComputeSingularityFactors() has been called.*/
std::complex<double> MGF::GetSingularityFactor(int component)
{
	if (!singularity_factors_computed)
	{
		throw std::logic_error("[ERROR] MGF::GetSingularityFactor(): Must call MGF::ComputeSingularityFactors() before requesting the singularity factors.");
	}

	if (component >= 0)
		return F[component];
	else
		return F_phi;
}


// ==================================================================================
// Computational drivers
// ==================================================================================

/*! \brief Driver to compute the MGF using numerical integration.*/
template<std::size_t N>
void MGF::ComputeMGF_Integration(double rho, double z, double zp, std::array<std::complex<double>, N> &G)
{

	std::fill(G.begin(), G.end(), 0.0);

	// If the wavelength in the source layer is extremely long, use only quasistatic terms
	if (s.extract_quasistatic)
	{
		bool quasistatic_only = UseQuasistaticOnly(rho, z, zp);
		if (quasistatic_only)
			return;
	}
	
	smgf.SetSourcePoint(zp);
	smgf.SetObservationPoint(z);

	if (s.components[0])
	{
		G[0] = IntegrateSpectralFarField(smgf, rho, 0, 0, s.switching_point);
		G[0] += IntegrateSpectralNearField(smgf, rho, 0, 0, s.switching_point);
	}
	
	if (s.components[1])
	{
		G[1] = IntegrateSpectralFarField(smgf, rho, 1, 1, s.switching_point);
		G[1] += IntegrateSpectralNearField(smgf, rho, 1, 1, s.switching_point);
	}

	if (s.components[2])
	{
		G[2] = IntegrateSpectralFarField(smgf, rho, 2, 1, s.switching_point);
		G[2] += IntegrateSpectralNearField(smgf, rho, 2, 1, s.switching_point);
	}

	if (s.components[3])
	{
		G[3] = IntegrateSpectralFarField(smgf, rho, 3, 0, s.switching_point);
		G[3] += IntegrateSpectralNearField(smgf, rho, 3, 0, s.switching_point);
	}

	if (s.components[4])
	{
		G[4] = IntegrateSpectralFarField(smgf, rho, 4, 0, s.switching_point);
		G[4] += IntegrateSpectralNearField(smgf, rho, 4, 0, s.switching_point);
	}
	
	return;

}


/*! \brief Driver to compute the curl MGF components using numerical integration.*/
template<std::size_t N>
void MGF::ComputeCurlMGF_Integration(double rho, double z, double zp, std::array<std::complex<double>, N> &G)
{

	std::fill(G.begin(), G.end(), 0.0);

	// If the wavelength in the source layer is extremely long, use only quasistatic terms
	if (s.extract_quasistatic)
	{
		bool quasistatic_only = UseQuasistaticOnly(rho, z, zp);		
		if (quasistatic_only)
			return;
	}
	
	smgf.SetSourcePoint(zp);
	smgf.SetObservationPoint(z);

	if (s.components_curl[0])
	{
		G[0] = IntegrateSpectralFarField(smgf, rho, 0, 2, s.switching_point, true);
		G[0] += IntegrateSpectralNearField(smgf, rho, 0, 2, s.switching_point, true);
	}
	
	if (s.components_curl[1])
	{
		G[1] = IntegrateSpectralFarField(smgf, rho, 1, 0, s.switching_point, true);
		G[1] += IntegrateSpectralNearField(smgf, rho, 1, 0, s.switching_point, true);
	}

	if (s.components_curl[2])
	{
		G[2] = IntegrateSpectralFarField(smgf, rho, 2, 1, s.switching_point, true);
		G[2] += IntegrateSpectralNearField(smgf, rho, 2, 1, s.switching_point, true);
	}

	if (s.components_curl[3])
	{
		G[3] = IntegrateSpectralFarField(smgf, rho, 3, 1, s.switching_point, true);
		G[3] += IntegrateSpectralNearField(smgf, rho, 3, 1, s.switching_point, true);
	}

	return;

}


/*! \brief Driver to compute the MGF using DCIM.*/
template<std::size_t N>
void MGF::ComputeMGF_DCIM(double rho, double z, double zp, std::array<std::complex<double>, N> &G)
{

	std::fill(G.begin(), G.end(), 0.0);

	// If the wavelength in the source layer is extremely long, use only quasistatic terms
	if (s.extract_quasistatic)
	{
		bool quasistatic_only = UseQuasistaticOnly(rho, z, zp);
		if (quasistatic_only)
			return;
	}
	
	std::vector<DCIM_images> im = dcim.GetImages(z, zp, false);
	double rho_sq = std::pow(rho, 2);

	// Traverse MGF components
	for (int ii = 0; ii < s.components.size(); ii++)
	{

		if (!s.components[ii])
			continue;

		std::complex<double> k = im[ii].k;

		// Traverse DCIM images
		for (int jj = 0; jj < im[ii].a.size(); jj++)
		{
			std::complex<double> a = im[ii].a[jj];
			std::complex<double> Rc = std::sqrt(rho_sq - std::pow(im[ii].alpha[jj], 2));

			if (ii == 0 || ii == 3 || ii == 4)
				G[ii] += a*std::exp(-J*k*Rc)/Rc;
			else
				G[ii] += rho*a*((1.0 + J*k*Rc)*std::exp(-J*k*Rc)/(std::pow(Rc, 3)));
		}

		G[ii] *= J/2.0/M_PI;

	}

	return;

}


/*! \brief Driver to compute the curl of the MGF using DCIM.*/
template<std::size_t N>
void MGF::ComputeCurlMGF_DCIM(double rho, double z, double zp, std::array<std::complex<double>, N> &G)
{

	std::fill(G.begin(), G.end(), 0.0);

	// If the wavelength in the source layer is extremely long, use only quasistatic terms
	if (s.extract_quasistatic)
	{
		bool quasistatic_only = UseQuasistaticOnly(rho, z, zp);
		if (quasistatic_only)
			return;
	}
	
	std::vector<DCIM_images> im = dcim.GetImages(z, zp, true);
	double rho_sq = std::pow(rho, 2);

	// Traverse MGF components
	for (int ii = 0; ii < s.components_curl.size(); ii++)
	{

		if (!s.components_curl[ii])
			continue;

		std::complex<double> k = im[ii].k;

		// Traverse DCIM images
		for (int jj = 0; jj < im[ii].a.size(); jj++)
		{
			std::complex<double> a = im[ii].a[jj];
			std::complex<double> Rc = std::sqrt(rho_sq - std::pow(im[ii].alpha[jj], 2));

			if (ii == 1)
				G[ii] += a*std::exp(-J*k*Rc)/Rc;
			else if (ii == 2 || ii == 3)
				G[ii] += rho*a*((1.0 + J*k*Rc)*std::exp(-J*k*Rc)/(std::pow(Rc, 3)));
			else
				G[ii] += a*std::pow(rho, 2)*((3.0 + 3.0*J*k*Rc + std::pow(J*k*Rc, 2))*std::exp(-J*k*Rc)/(std::pow(Rc, 5)));		
		}

		G[ii] *= J/2.0/M_PI;

	}

	return;

}


/*! \brief Driver to compute the MGF using a pre-computed interpolation table.*/
template<std::size_t N>
void MGF::ComputeMGF_Interpolation(double rho, double z, double zp, std::array<std::complex<double>, N> &G, std::vector<std::vector<table_entry<N>>> &table, std::vector<bool> &components)
{

	std::fill(G.begin(), G.end(), 0.0);

	// Extract interpolation abscissae
	int row = GetRow(z, zp);
	std::vector<int> cols = GetColumns(rho);

	// Interpolate the MGF
	for (int qq = 0; qq < cols.size(); qq++)
	{

		// Compute the Lagrange polynomials
		double L = 1.0;
		for (int ii = 0; ii < cols.size(); ii++)
		{
			if (ii == qq)
				continue;
			L *= (rho - lm.rho_nodes[cols[ii]])/(lm.rho_nodes[cols[qq]] - lm.rho_nodes[cols[ii]]);
		}

		for (int ii = 0; ii < G.size(); ii++)
			if (components[ii])
				G[ii] += table[row][cols[qq]].K[ii]*L;
		
	}

	return;

}


/*! \brief Function to compute the singularity factors for the quasistatic MGF.*/
void MGF::ComputeSingularityFactors(double x_diff, double y_diff, double z, double zp)
{

	double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

	std::complex<double> zeta = std::atan2(y_diff, x_diff);
	cos_term = std::cos(zeta);
	sin_term = std::sin(zeta);

	qmgf.ComputeSingularityFactors(z, zp, rho);
	ComputeSingularityFactors();

	return;
	
}


/*! \brief Function to compute the singularity factors for the quasistatic MGF, assuming that ComputeMGF() has already been called.*/
void MGF::ComputeSingularityFactors()
{

	std::fill(F.begin(), F.end(), 0.0);

	F[0] = qmgf.GetSingularityFactor(0);
	F[1] = 0.0;
	F[2] = qmgf.GetSingularityFactor(1)*cos_term;
	F[3] = 0.0;
	F[4] = qmgf.GetSingularityFactor(0);
	F[5] = qmgf.GetSingularityFactor(1)*sin_term;
	F[6] = qmgf.GetSingularityFactor(2)*cos_term;
	F[7] = qmgf.GetSingularityFactor(2)*sin_term;
	F[8] = qmgf.GetSingularityFactor(3);
	
	F_phi = qmgf.GetSingularityFactor(4);

	singularity_factors_computed = true;
	
	return;

}


/*! \brief Function to compute the homogeneous term factors for the MGF.*/
void MGF::ComputeHomogeneousFactors()
{

	std::fill(F.begin(), F.end(), 0.0);

	std::complex<double> epsi, mui;

	if (i == -1)
	{
		epsi = lm.eps_top;
		mui = lm.mu_top;
	}
	else if (i == lm.layers.size())
	{
		epsi = lm.eps_bot;
		mui = lm.mu_bot;
	}
	else
	{
		epsi = lm.eps[i];
		mui = lm.mu[i];
	}
	
	F[0] = mui/mu0/(4.0*M_PI);
	F[4] = mui/mu0/(4.0*M_PI);
	F[8] = mui/mu0/(4.0*M_PI);
	F_phi = eps0/epsi/(4.0*M_PI);

	singularity_factors_computed = true;
	
	return;

}


// ==================================================================================
// Computational helpers
// ==================================================================================

/*! \brief Check if distance is small enough that quasistatic terms alone are sufficient.*/
bool MGF::UseQuasistaticOnly(double rho, double z, double zp)
{

	bool quasistatic_only = false;
	double r = std::sqrt(std::pow(rho, 2) + std::pow((z-zp), 2));
	std::complex<double> g0 = std::exp(-J*lm.k[i]*r);
	// if (std::abs(1.0 - std::abs(std::real(g0))) < s.tol_qse && std::abs(std::imag(g0)) < s.tol_qse)
	if (std::abs(std::real(lm.k[i]*r)) < s.tol_qse)
		quasistatic_only = true;

	// temp
	// return true;
	// if (std::abs(r) < 5.0e-6)
	// 	quasistatic_only = true;
	
	return quasistatic_only;

}


/*! \brief Generate an interpolation table for the spatial MGF or its curl components.*/
template<std::size_t N>
void MGF::TabulateMGF(std::vector<std::vector<table_entry<N>>> &table, bool curl)
{

	if (lm.z_nodes.size() < 1)
	{
		std::cout << "[WARNING] MGF::TabulateMGF(): No z-nodes have been defined, so no samples were tabulated." << std::endl;
		return;
	}

	if (lm.rho_nodes.size() < 1)
	{
		std::cout << "[WARNING] MGF::TabulateMGF(): No rho-nodes have been defined, so no samples were tabulated." << std::endl;
		return;
	}


	// ====== Reset table ======

    table.clear();


	// ====== Generate index maps ======

	GenerateTableMaps();
	

	// ====== Generate table ======

	table.reserve(z_to_idx.size());

	// Traverse source layers
	for (int ii = 0; ii < lm.layers.size(); ii++)
	{
		// Traverse observer layers
		for (int mm = 0; mm < lm.layers.size(); mm++)
		{
			smgf.SetLayers(ii, mm);
			i = ii;
			m = mm;
			
			// Traverse source z-nodes
			for (int ss = 0; ss < lm.z_nodes[ii].size(); ss++)
			{
				// Traverse observer z-nodes
				for (int tt = 0; tt < lm.z_nodes[mm].size(); tt++)
				{

					double zp = lm.z_nodes[ii][ss];
					double z = lm.z_nodes[mm][tt];					

					// Generate all entries for this row
					table.push_back(std::vector<table_entry<N>> (lm.rho_nodes.size()));

					for (int qq = 0; qq < lm.rho_nodes.size(); qq++)
					{
						double rho = lm.rho_nodes[qq];

						if (!curl)
						{
							if (s.sampling_method == MGF_INTEGRATE)
								ComputeMGF_Integration(rho, z, zp, table.back()[qq].K);
							else if (s.sampling_method == MGF_DCIM)
								ComputeMGF_DCIM(rho, z, zp, table.back()[qq].K);
						}
						else
						{
							if (s.sampling_method == MGF_INTEGRATE)
								ComputeCurlMGF_Integration(rho, z, zp, table.back()[qq].K);
							else if (s.sampling_method == MGF_DCIM)
								ComputeCurlMGF_DCIM(rho, z, zp, table.back()[qq].K);
						}
					}
				}
			}
		}
	}

	return;

}


/*! \brief Generate maps for interpolation table access.*/
void MGF::GenerateTableMaps()
{

	z_to_idx.clear();
	rho_to_idx.clear();
	idxpair_to_row.clear();
	

	// ====== Generate maps between z-nodes and their index in the stackup ======

	// The int pair <ii, zz> is mapped to a unique int using Szudzik's function.
	// Reference: https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
	for (int ii = 0; ii < lm.z_nodes.size(); ii++)
	{
		for (int zz = 0; zz < lm.z_nodes[ii].size(); zz++)
		{
			int idx = ii >= zz ? ii*ii + ii + zz : ii + zz*zz;
			z_to_idx.insert(std::make_pair(lm.z_nodes[ii][zz], idx));
		}
	}

	int map_size = z_to_idx.size();


	// ====== Generate a map for rho-nodes ======

	for (int ii = 0; ii < lm.rho_nodes.size(); ii++)
		rho_to_idx.insert(std::make_pair(lm.rho_nodes[ii], ii));


	// ====== Generate a map for table entries ======
	
	int idx_row = 0;
	
	// Traverse source layers
	for (int ii = 0; ii < lm.layers.size(); ii++)
	{
		// Traverse observer layers
		for (int mm = 0; mm < lm.layers.size(); mm++)
		{
			// Traverse source z-nodes
			for (int ss = 0; ss < lm.z_nodes[ii].size(); ss++)
			{
				// Traverse observer z-nodes
				for (int tt = 0; tt < lm.z_nodes[mm].size(); tt++)
				{

					double zp = lm.z_nodes[ii][ss];
					double z = lm.z_nodes[mm][tt];

					int idx_z = z_to_idx[z];
					int idx_zp = z_to_idx[zp];

					idxpair_to_row.insert(std::make_pair(std::make_pair(idx_z, idx_zp), idx_row));

					idx_row++;

				}
			}
		}
	}

	int map_size_new = z_to_idx.size();

	if (map_size_new != map_size)
		std::cout << "[WARNING] MGF::GenerateTableMaps(): Table indexing may have been corrupted due to numerical issues." << std::endl;	

	
	return;

}


/*! \brief Function to retrieve the nearest tabulated row for a given z-zp pair.*/
int MGF::GetRow(double z, double zp)
{

	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF:GetRow(): Table has not been generated. Call MGF::Initialize() first.");
	}


	// ====== Locate the index of the z-node nearest to z ======

	std::map<double, int>::iterator z_below = z_to_idx.lower_bound(z);
	std::map<double, int>::iterator z_above = z_to_idx.upper_bound(z);

	int idx_z;
	double z_found;

	if (std::abs(z - z_below->first) < 1.0e-15)
	{
		idx_z = z_below->second;
		z_found = z_below->first;
	}
	else
	{
		z_below--;
	
		if (std::abs(z - z_below->first) < std::abs(z - z_above->first))
		{
			idx_z = z_below->second;
			z_found = z_below->first;
		}
		else
		{
			idx_z = z_above->second;
			z_found = z_above->first;
		}
	}


	// ====== Locate the index of the z-node nearest to zp ======

	std::map<double, int>::iterator zp_below = z_to_idx.lower_bound(zp);
	std::map<double, int>::iterator zp_above = z_to_idx.upper_bound(zp);

	int idx_zp;
	double zp_found;

	if (std::abs(zp - zp_below->first) < 1.0e-15)
	{
		idx_zp = zp_below->second;
		zp_found = zp_below->first;
	}
	else
	{
		zp_below--;
	
		if (std::abs(zp - zp_below->first) < std::abs(zp - zp_above->first))
		{
			idx_zp = zp_below->second;
			zp_found = zp_below->first;
		}
		else
		{
			idx_zp = zp_above->second;
			zp_found = zp_above->first;
		}
	}
	

	// ====== Retrive the images at the above index pair ======

	int idx_row = idxpair_to_row[std::make_pair(idx_z, idx_zp)];

	return idx_row;

}


/*! \brief Function to retrieve the stencil points along rows for a given rho value.*/
std::vector<int> MGF::GetColumns(double rho)
{

	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF:GetColumns(): Table has not been generated. Call MGF::Initialize() first.");
	}

	std::vector<int> rho_stencil (s.order+1);
	

	// ====== Locate the index of the rho-node nearest to rho ======

	std::map<double, int>::iterator rho_below = rho_to_idx.lower_bound(rho);
	std::map<double, int>::iterator rho_above = rho_to_idx.upper_bound(rho);
	std::map<double, int>::iterator rho_it;

	if (std::abs(rho - rho_below->first) < 1.0e-15)
		rho_it = rho_below;
	else
	{
		rho_below--;
	
		if (std::abs(rho - rho_below->first) < std::abs(rho - rho_above->first))
			rho_it = rho_below;
		else
			rho_it = rho_above;
	}

	int idx_rho = std::distance(rho_to_idx.begin(), rho_it);

	int idx_start = std::max((idx_rho - s.order/2), 0);
	if (idx_start + s.order > rho_to_idx.size() - 1)
		idx_start = rho_to_idx.size() - s.order - 1;

	int start_dist = idx_rho - idx_start;
	std::advance(rho_it, -start_dist);
	
	for (int ii = idx_start, jj = 0; ii <= idx_start + s.order; ii++, jj++)
	{
		rho_stencil[jj] = rho_it->second;
		rho_it++;
	}

	// ====== Debugging ======

	// std::cout << rho_it->first << ", " << idx_rho << ", " << idx_start << ", " << start_dist << std::endl;

	// std::cout << rho << ", [ " << std::flush;
	// for (int ii = 0; ii < rho_stencil.size(); ii++)
	// 	std::cout << lm.rho_nodes[rho_stencil[ii]] << " " << std::flush;
	// std::cout << "]" << std::endl;

	// std::ofstream print_rho ("./test/benchmarks/MGF_layers/test_rho.debug", std::ios_base::app);
	// print_rho << rho << ", [ " << std::flush;
	// for (int ii = 0; ii < rho_stencil.size(); ii++)
	// 	print_rho << lm.rho_nodes[rho_stencil[ii]] << " " << std::flush;
	// print_rho << "]" << std::endl;

	return rho_stencil;

}


/*! \brief Load an MGF interpolation table from disk. Returns 0 if successful.*/
template<std::size_t N>
int MGF::LoadTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename)
{

	std::ifstream fin;
	fin.open(filename, std::ios::binary);

	if (!fin)
	{
		std::cout << "\n[WARNING] MGF::LoadTable(): Failed to open the interpolation table file. Skipping import and generating a new table instead." << std::endl;
		return -1;
	}

	// Read in total number of z_nodes
	int N_z_nodes;
	fin.read(reinterpret_cast<char *> (&N_z_nodes), sizeof(N_z_nodes));

	// Read in z-nodes and the layers to which they belong
	lm.ClearNodes_z();
	for (int ii = 0; ii < N_z_nodes; ii++)
	{
		struct X
		{
			int layer;
			double z;
		} x;

		fin.read(reinterpret_cast<char *> (&x), sizeof(x));

		std::vector<double> z_node (1);
		z_node[0] = x.z;
		lm.InsertNodes_z(z_node, x.layer);
	}

	// Read in total number of rho_nodes
	int N_rho_nodes;
	fin.read(reinterpret_cast<char *> (&N_rho_nodes), sizeof(N_rho_nodes));

	// Read in rho-nodes
	std::vector<double> rho_nodes (N_rho_nodes);
	for (int ii = 0; ii < N_rho_nodes; ii++)
	{
		double rho;
		fin.read(reinterpret_cast<char *> (&rho), sizeof(rho));
		rho_nodes[ii] = rho;
	}
	lm.ClearNodes_rho();
	lm.InsertNodes_rho(rho_nodes);

	// Read in the MGF table
	table.clear();
	table.resize(N_z_nodes*N_z_nodes);
	for (int ii = 0; ii < N_z_nodes*N_z_nodes; ii++)
	{
		table[ii].resize(N_rho_nodes);
		
		for (int jj = 0; jj < N_rho_nodes; jj++)
		{
			for (int kk = 0; kk < table[ii][jj].K.size(); kk++)
				fin.read(reinterpret_cast<char *> (&table[ii][jj].K[kk]), sizeof(table[ii][jj].K[kk]));
		}
	}

	// Generate index maps
	GenerateTableMaps();
	
	return 0;

}


/*! \brief Export an MGF interpolation table to disk. Returns 0 if successful.*/
template<std::size_t N>
int MGF::ExportTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename)
{

	std::ofstream fout;
	fout.open(filename, std::ios::binary);

	if (!fout)
	{
		std::cout << "\n[WARNING] MGF::ExportTable(): Failed to create or open the interpolation table file. Skipping export." << std::endl;
		return -1;
	}

	// Write out header to record the total number of z_nodes
	int N_z_nodes = 0;
	for (int ii = 0; ii < lm.z_nodes.size(); ii++)
		N_z_nodes += lm.z_nodes[ii].size();
	fout.write(reinterpret_cast<char *> (&N_z_nodes), sizeof(N_z_nodes));

	// Write out z-nodes and the layers to which they belong
	for (int ii = 0; ii < lm.z_nodes.size(); ii++)
	{
		for (int jj = 0; jj < lm.z_nodes[ii].size(); jj++)
		{
			struct X
			{
				int layer;
				double z;
			} x;

			x.layer = ii;
			x.z = lm.z_nodes[ii][jj];
			fout.write(reinterpret_cast<char *> (&x), sizeof(x));
		}
	}

	// Write out header to record the total number of rho_nodes
	int N_rho_nodes = lm.rho_nodes.size();
	fout.write(reinterpret_cast<char *> (&N_rho_nodes), sizeof(N_rho_nodes));

	// Write out rho-nodes
	for (int ii = 0; ii < lm.rho_nodes.size(); ii++)
	{
		double rho = lm.rho_nodes[ii];
		fout.write(reinterpret_cast<char *> (&rho), sizeof(rho));
	}
	
	// Write out the MGF table
	for (int ii = 0; ii < table.size(); ii++)
	{
		for (int jj = 0; jj < table[ii].size(); jj++)
		{
			for (int kk = 0; kk < table[ii][jj].K.size(); kk++)
				fout.write(reinterpret_cast<char *> (&table[ii][jj].K[kk]), sizeof(table[ii][jj].K[kk]));
		}
	}
	

	return 0;

}


