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
#include <numeric>

#include "MGF.hpp"
#include "layers.hpp"
#include "spectral_MGF.hpp"
#include "quasistatic_MGF.hpp"
#include "sommerfeld_integrals.hpp"
#include "DCIM.hpp"
#include "constants.hpp"
#include "progress_bar.hpp"

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

/*! \brief Consider one layer first*/
void MGF::adaptiveInterpolation(MGF &mgf)
{
    updateRhonodes(mgf);
    //updateZnodes(mgf);

    return;
}

/*! \brief Consider one layer first*/
void MGF::updateZnodes(MGF &mgf)
{

    std::cout << "========================= Update the z nodes =========================" << std::endl;

    for (int layer = 0; layer < mgf.lm.layers.size(); layer++)
    {
        // Ensure the vector has at least 3 elements
        if (mgf.lm.z_nodes[layer].size() < 2) {
            std::cerr << "Vector needs at least 2 elements." << std::endl;
        }

        std::vector<double> z_test_nodes;

        for (int ii = 0; ii < mgf.lm.z_nodes[layer].size() - 1; ii++)
        {
            double z_test = (mgf.lm.z_nodes[layer][ii+1] + mgf.lm.z_nodes[layer][ii])/2;
            z_test_nodes.push_back(z_test);
        }

        // compute the MGF
        // choose rho that will more likely to lead to strong singularity
        double rho_test = mgf.lm.rho_nodes[1];

        for (int ii = 0; ii < z_test_nodes.size(); ii++)
        {
            double z_test = z_test_nodes[ii];
            double z_src = z_test_nodes[ii];

            double z_spacing = mgf.lm.z_nodes[layer][ii + 1] - mgf.lm.z_nodes[layer][ii];
            int level = 1;

            int i = lm.FindLayer(z_test);
            int m = lm.FindLayer(z_src);
            mgf.smgf.SetLayers(i, m);
            bool is_Midpoint_Correct = isMidpointCorrect(rho_test, z_src, z_test, mgf.s.adaptive_threshold, false);

            if (!is_Midpoint_Correct)
            {
                if(level == 1)
                    addZTable(mgf, z_test, layer);
                level++;
                std::vector<double> z_tests;
                test_addZTable_recursive(mgf, layer, rho_test, z_spacing, z_test, z_src, level, z_tests);
                processZTests(mgf, layer, rho_test, z_spacing, z_test, z_src, level, z_tests);
            }
        }

    }

    return;
}

/*! \brief Consider one layer first*/
void MGF::updateRhonodes(MGF &mgf)
{
    std::cout << "===========================Update the rho nodes===========================" << std::endl;

    std::vector<double> rho_test_nodes;

	std::cout << "===========================Initialize the test rho nodes===========================" << std::endl;
    for (int ii = 0; ii < mgf.lm.rho_nodes.size() - 1; ii++)
    {
        double rho_test = (mgf.lm.rho_nodes[ii] + mgf.lm.rho_nodes[ii+1])/2;
        rho_test_nodes.push_back(rho_test);
    }


    // compute the MGF
    // choose rho that will more likely to lead to strong singularity

    double z_test = mgf.lm.z_nodes[0][0];
    double z_src = mgf.lm.z_nodes[0][0];

	std::cout << "===========================Iteration on the test rho nodes===========================" << std::endl;
    for (int ii = 0; ii < rho_test_nodes.size(); ii++) {

        double rho_test = rho_test_nodes[ii];
        double rho_spacing = lm.rho_nodes[ii + 1] - lm.rho_nodes[ii];
        int level = 1;

        // verify if this point will be added based on the relative error rate
        int i = lm.FindLayer(z_test);
        int m = lm.FindLayer(z_src);
        mgf.smgf.SetLayers(i, m);
        mgf.i = i;
        mgf.m = m;


    	std::cout << "===========================Initialization before is_Midpoint_Correct===========================" << std::endl;
        bool is_Midpoint_Correct = isMidpointCorrect(rho_test, z_test, z_src, mgf.s.adaptive_threshold, true);

        if (!is_Midpoint_Correct)
        {
            if(level == 1)
                addRhoTable(mgf, rho_test);
            level++;
            std::vector<double> rho_test_nodes_new;
            test_addRhoTable_recursive(mgf, rho_test, rho_spacing, z_test, z_src, level, rho_test_nodes_new);
            processRhoTests(mgf, rho_spacing, z_test, z_src, level, rho_test_nodes_new);
        }
    }


    return;
}

void MGF::processRhoTests(MGF& mgf, double rho_spacing, double z_test, double z_src, int level, std::vector<double>& rho_tests)
{
    if (level >= 6)
        return;
    for (double test_rho : rho_tests) {
        std::vector<double> new_rho_tests;
        test_addRhoTable_recursive(mgf, test_rho, rho_spacing, z_test, z_src, level + 1, new_rho_tests);

        // Process new test points if generated.
        if (!new_rho_tests.empty()) {
            processRhoTests(mgf, rho_spacing, z_test, z_src, level + 1, new_rho_tests);
        }
    }
}

void MGF::processZTests(MGF& mgf, int layer_idx, double rho_test, double z_spacing, double z_test, double z_src, int level, std::vector<double>& z_tests)
{
    if (level >= 3)
        return;
    for (double test_z : z_tests) {
        std::vector<double> new_z_tests;
        test_addZTable_recursive(mgf, layer_idx, rho_test, z_spacing, z_test, z_src, level + 1, new_z_tests);

        // Process new test points if generated.
        if (!new_z_tests.empty()) {
            processZTests(mgf, layer_idx, rho_test, z_spacing, test_z, test_z, level + 1, new_z_tests);
        }
    }
}


bool MGF::isMidpointCorrect(double rho, double z_src, double z_test, double adaptive_threshold, bool test_rho)
{
    std::array<std::complex<double>, 5> _G_integ;
    std::array<std::complex<double>, 5> _G_interp;

    std::fill(_G_integ.begin(), _G_integ.end(), 0.0);
    std::fill(_G_interp.begin(), _G_interp.end(), 0.0);

	std::cout << "===========================ComputeMGF_Integration===========================" << std::endl;

    ComputeMGF_Integration(rho, z_test, z_src, _G_integ);

	std::cout << "===========================ComputeMGF_Interpolation_withZ===========================" << std::endl;

    ComputeMGF_Interpolation_withZ(rho, z_test, z_src, _G_interp, MGF_table, s.components);

    std::vector<double> rmse;
    for (int i = 0; i < 5; i++)
    {
        //std::cout << i << ": G_integration = " << _G_integ[i] << ", G_interplation = " << _G_interp[i] << std::endl;

        if (std::abs(_G_integ[i]) != 0)
            rmse.push_back(std::abs(_G_integ[i] - _G_interp[i]) / std::abs(_G_integ[i]));
    }

    double rmse_mean = std::accumulate(rmse.begin(),rmse.end(),0.0) / rmse.size();

    if (test_rho)
        std::cout << "The RMSE for interpolation point at rho = " << rho << " is " << rmse_mean << std::endl;
    else
        std::cout << "The RMSE for interpolation point at z = " << z_test << " is " << rmse_mean << std::endl;
    if (rmse_mean > adaptive_threshold)
    {
        return false;
    }

    return true;
}
void MGF::addRhoTable(MGF &mgf, double rho_test)
{
    int rho_size = mgf.lm.rho_nodes.size();
    // Find the proper position to insert the value so that the vector remains sorted
    auto position = std::lower_bound(lm.rho_nodes.begin(), lm.rho_nodes.end(), rho_test);

    // Calculate the index for the new element
    std::vector<int>::difference_type index = position - lm.rho_nodes.begin();

    // Insert the value
    mgf.lm.rho_nodes.insert(position, rho_test);

    // Update the interpolation table
    AppendMGFTable_rho(mgf.MGF_table, index, rho_size + 1);

    std::cout << "-->Add this point to the interpolation table!" << std::endl;
}

void MGF::addZTable(MGF &mgf, double z_test, int layer)
{
    int z_size = mgf.lm.z_nodes[layer].size();
    // Find the proper position to insert the value so that the vector remains sorted
    auto position = std::lower_bound(mgf.lm.z_nodes[layer].begin(), mgf.lm.z_nodes[layer].end(), z_test);

    // Calculate the index for the new element
    std::vector<int>::difference_type index = position - mgf.lm.z_nodes[layer].begin();

    // Insert the value
    mgf.lm.z_nodes[layer].insert(position, z_test);

    // Update the interpolation table
    AppendMGFTable_z(mgf.MGF_table, layer, index, z_size + 1);

    std::cout << "-->Add this point to the interpolation table!" << std::endl;
}

void MGF::test_addRhoTable_recursive(MGF &mgf, double rho_test, double rho_spacing, double z_test, double z_src, int level, std::vector<double> &rho_tests_l2)
{
    double spacing_l1 = rho_spacing / std::pow(2, level);

    std::vector<double> rho_tests_l1 {rho_test - spacing_l1, rho_test + spacing_l1};

    for (int jj = 0; jj < 2; jj++)
    {
        int i = lm.FindLayer(z_test);
        int m = lm.FindLayer(z_src);
        mgf.smgf.SetLayers(i, m);

        bool add_point_to_table = isMidpointCorrect(rho_tests_l1[jj], z_test, z_src, mgf.s.adaptive_threshold, true);
        if (!add_point_to_table)
        {
            addRhoTable(mgf, rho_tests_l1[jj]);
            rho_tests_l2.push_back(rho_tests_l1[jj]);
        }
    }
}

void MGF::test_addZTable_recursive(MGF &mgf, int layer_idx, double rho_test, double z_spacing, double z_test, double z_src, int level, std::vector<double> &z_tests_l2)
{
    double spacing_l1 = z_spacing / std::pow(2, level);

    std::vector<double> z_tests_l1 {z_test - spacing_l1, z_test + spacing_l1};

    for (int jj = 0; jj < 2; jj++)
    {
        int i = lm.FindLayer(z_tests_l1[jj]);
        int m = lm.FindLayer(z_tests_l1[jj]);
        mgf.smgf.SetLayers(i, m);
        bool is_Midpoint_Correct = isMidpointCorrect(rho_test, z_tests_l1[jj], z_tests_l1[jj], mgf.s.adaptive_threshold, false);

        if (!is_Midpoint_Correct)
        {
            addZTable(mgf, z_tests_l1[jj], layer_idx);
            z_tests_l2.push_back(z_tests_l1[jj]);
        }
    }
}


void MGF::plotRhoNodes(MGF mgf, std::vector<double> z_gridpoints)
{
    double z_test = z_gridpoints[0];
    double z_src = z_gridpoints[0];

    std::ofstream outputFile("../Testing/MGF.txt");
    for (int ii = 0; ii < mgf.lm.rho_nodes.size(); ii++)
    {
        outputFile << mgf.lm.rho_nodes[ii] << ", ";
        std::array<std::complex<double>, 5> _G_integ;
        std::fill(_G_integ.begin(), _G_integ.end(), 0.0);
        ComputeMGF_Integration(mgf.lm.rho_nodes[ii], z_test, z_src, _G_integ);

        if(outputFile.is_open()){
            for (int jj = 0; jj < _G_integ.size(); jj++)
            {
                outputFile << _G_integ[jj];
                if(jj != _G_integ.size() - 1)
                    outputFile << ", ";
                if(jj == _G_integ.size() - 1)
                    outputFile << "\n";
            }
        }
        else
        {
            std::cout << "Can't open the file!" << std::endl;
        }

    }
    outputFile.close();



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
        if(s.interpolate_z)
		    ComputeMGF_Interpolation_withZ(rho, z, zp, _G, MGF_table, s.components);
        else
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

/*! \brief Driver to compute the MGF using a pre-computed interpolation table.*/
template<std::size_t N>
void MGF::ComputeMGF_Interpolation_withZ(double rho, double z, double zp, std::array<std::complex<double>, N> &G, std::vector<std::vector<table_entry<N>>> &table, std::vector<bool> &components)
{

    std::fill(G.begin(), G.end(), 0.0);

    // Extract interpolation abscissae (index and position)
    std::vector<int> z_idx_stencil;
    std::vector<int> zp_idx_stencil;

    std::vector<double> z_stencil;
    std::vector<double> zp_stencil;

    GetStencil_z(z, z_idx_stencil, z_stencil);
    GetStencil_z(zp, zp_idx_stencil, zp_stencil);

    // Get interpolation points for rho
    std::vector<int> cols = GetColumns(rho);

    // Interpolate the MGF
    for (int ii = 0; ii < cols.size(); ii++)
    {
        for (int jj = 0; jj < z_stencil.size(); jj++)
        {
            for (int kk = 0; kk < zp_stencil.size(); kk++)
            {
                // Compute the Lagrange polynomials
                double L = 1.0;

                // ====== Retrive the images at the above index pair ======

                int idx_row = idxpair_to_row[std::make_pair(z_idx_stencil[jj], zp_idx_stencil[kk])];

                for (int tt = 0; tt < cols.size(); tt++)
                {
                    if (tt == ii)
                        continue;
                    L *= (rho - lm.rho_nodes[cols[tt]])/(lm.rho_nodes[cols[ii]] - lm.rho_nodes[cols[tt]]);
                }

                for (int tt = 0; tt < z_stencil.size(); tt++)
                {
                    if (tt == jj)
                        continue;
                    L *= (z - z_stencil[tt])/(z_stencil[jj] - z_stencil[tt]);
                }

                for (int tt = 0; tt < zp_stencil.size(); tt++)
                {
                    if (tt == kk)
                        continue;
                    L *= (zp - zp_stencil[tt])/(zp_stencil[kk] - zp_stencil[tt]);
                }


                for (int tt = 0; tt < G.size(); tt++)
                    if (components[tt])
                        G[tt] += table[idx_row][cols[ii]].K[tt]*L;
            }
        }

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

	// Precompute row‑block offsets for each (ii,mm) pair
	std::vector<int> layerOffsets(lm.layers.size() * lm.layers.size());
	int offset = 0;
	for (int ii = 0; ii < lm.layers.size(); ++ii) {
		for (int mm = 0; mm < lm.layers.size(); ++mm) {
			layerOffsets[ii * lm.layers.size() + mm] = offset;
			offset +=
				static_cast<int>(lm.z_nodes[ii].size()) *
				static_cast<int>(lm.z_nodes[mm].size());
		}
	}
	const int totalRows = offset;

	// Pre‑allocate the table: totalRows × lm.rho_nodes.size()
	table.assign(totalRows, std::vector<table_entry<N>>(lm.rho_nodes.size()));


	// ====== Generate index maps ======

	GenerateTableMaps();

	// ====== Generate table ======

	MGF mgfLocal = *this;
	// loop‐indices and temporaries for the parallel region
	int ii, mm, ss, tt, qq, rowIdx, base;
	double z, zp, rho;


#pragma omp parallel for collapse(2) firstprivate(mgfLocal) schedule(dynamic) \
    default(none) \
    shared(table, lm, s, layerOffsets, curl) \
    private(ii, mm, ss, tt, qq, rowIdx, rho, z, zp, base)

	// Traverse source layers
	for (ii = 0; ii < lm.layers.size(); ii++)
	{

		// Traverse observer layers
		for (mm = 0; mm < lm.layers.size(); mm++)
		{
			base = layerOffsets[ii*lm.layers.size() + mm];
			int zLenSrc = (int)lm.z_nodes[ii].size();
			int zLenObs = (int)lm.z_nodes[mm].size();

			// Traverse source z-nodes
			for (ss = 0; ss < zLenSrc; ss++)
			{
				// Traverse observer z-nodes
				for (tt = 0; tt < zLenObs; tt++)
				{

					// Compute the unique row index
					rowIdx = base + ss*zLenObs + tt;

					// Update thread‑local state
					mgfLocal.i = ii;
					mgfLocal.m = mm;
					mgfLocal.smgf.SetLayers(ii, mm);

					zp = lm.z_nodes[ii][ss];
					z = lm.z_nodes[mm][tt];


					for (qq = 0; qq < lm.rho_nodes.size(); qq++)
					{
						rho = lm.rho_nodes[qq];
						auto &Kref = table[rowIdx][qq].K;

						if (!curl)
						{
							if (s.sampling_method == MGF_INTEGRATE)
								mgfLocal.ComputeMGF_Integration(rho, z, zp, Kref);
							else if (s.sampling_method == MGF_DCIM)
								mgfLocal.ComputeMGF_DCIM(rho, z, zp, Kref);
						}
						else
						{
							if (s.sampling_method == MGF_INTEGRATE)
								mgfLocal.ComputeCurlMGF_Integration(rho, z, zp, Kref);
							else if (s.sampling_method == MGF_DCIM)
								mgfLocal.ComputeCurlMGF_DCIM(rho, z, zp, Kref);
						}
					}
				}
			}
		}
	}

	return;

}

template<std::size_t N>
void MGF::AppendMGFTable_z(std::vector<std::vector<table_entry<N>>> &table, int layer_idx, int z_idx, int z_new_size, bool curl)
{

    if (lm.z_nodes.size() < 1)
    {
        std::cout << "[WARNING] MGF::AppendMGFTable(): No z-nodes have been defined, so no samples were tabulated." << std::endl;
        return;
    }


    // ====== Generate index maps ======

    AddTableMaps_z(layer_idx, z_idx, z_new_size);


    // ====== Generate table ======

    // Traverse observer z-nodes
    for (int ii = 0; ii < lm.layers.size(); ii++)
    {
        smgf.SetLayers(layer_idx, ii);

        for (int tt = 0; tt < lm.z_nodes[ii].size(); tt++) {

            double zp = lm.z_nodes[layer_idx][z_idx];
            double z = lm.z_nodes[ii][tt];

            // Generate all entries for this row
            table.push_back(std::vector<table_entry<N>>(lm.rho_nodes.size()));

            for (int qq = 0; qq < lm.rho_nodes.size(); qq++) {
                double rho = lm.rho_nodes[qq];

                if (!curl) {
                    if (s.sampling_method == MGF_INTEGRATE)
                        ComputeMGF_Integration(rho, z, zp, table.back()[qq].K);
                    else if (s.sampling_method == MGF_DCIM)
                        ComputeMGF_DCIM(rho, z, zp, table.back()[qq].K);
                } else {
                    if (s.sampling_method == MGF_INTEGRATE)
                        ComputeCurlMGF_Integration(rho, z, zp, table.back()[qq].K);
                    else if (s.sampling_method == MGF_DCIM)
                        ComputeCurlMGF_DCIM(rho, z, zp, table.back()[qq].K);
                }
            }
        }
    }


    // Traverse source z-nodes
    for (int ii = 0; ii < lm.layers.size(); ii++)
    {
        smgf.SetLayers(ii, layer_idx);

        for (int tt = 0; tt < lm.z_nodes[ii].size(); tt++)
        {
            double zp = lm.z_nodes[ii][tt];

            if (zp == lm.z_nodes[layer_idx][z_idx])
                continue;

            double z = lm.z_nodes[layer_idx][z_idx];

            // Generate all entries for this row
            table.push_back(std::vector<table_entry<N>>(lm.rho_nodes.size()));

            for (int qq = 0; qq < lm.rho_nodes.size(); qq++) {
                double rho = lm.rho_nodes[qq];

                if (!curl) {
                    if (s.sampling_method == MGF_INTEGRATE)
                        ComputeMGF_Integration(rho, z, zp, table.back()[qq].K);
                    else if (s.sampling_method == MGF_DCIM)
                        ComputeMGF_DCIM(rho, z, zp, table.back()[qq].K);
                } else {
                    if (s.sampling_method == MGF_INTEGRATE)
                        ComputeCurlMGF_Integration(rho, z, zp, table.back()[qq].K);
                    else if (s.sampling_method == MGF_DCIM)
                        ComputeCurlMGF_DCIM(rho, z, zp, table.back()[qq].K);
                }
            }
        }
    }



    return;

}

template<std::size_t N>
void MGF::AppendMGFTable_rho(std::vector<std::vector<table_entry<N>>> &table, int rho_idx, int rho_new_size, bool curl)
{

    if (lm.rho_nodes.size() < 1)
    {
        std::cout << "[WARNING] MGF::AppendMGFTable(): No rho-nodes have been defined, so no samples were tabulated." << std::endl;
        return;
    }


    // ====== Generate index maps ======

    AddTableMaps_rho();

	// Grow every row to make room for one more rho‑column
    for (auto &row : table)
        row.resize(rho_new_size);

	// Pull out the fixed number of layers and rhos
    const int L = static_cast<int>(lm.layers.size());
    const int R = static_cast<int>(lm.rho_nodes.size());


    // Make a thread‑local copy of *this* (carries i,m, smgf, maps, settings…)
    MGF mgfLocal = *this;
	int ii, mm, ss, tt, idx_z, idx_zp, idx_row;
    int zLenSrc, zLenObs;
    double z, zp, rho;

    #pragma omp parallel for collapse(2) \
    firstprivate(mgfLocal, rho_idx, rho_new_size, curl, L, R) \
    schedule(dynamic) default(none) \
    shared(table) \
    private(ii, mm, ss, tt, idx_z, idx_zp, idx_row, zLenSrc, zLenObs, z, zp, rho)
	// Traverse source layers
    for (ii = 0; ii < L; ii++)
    {
    	// Traverse observer layers
        for (mm = 0; mm < L; mm++)
        {
        	// Cache these per (ii,mm)
			zLenSrc = static_cast<int>(lm.z_nodes[ii].size());
			zLenObs = static_cast<int>(lm.z_nodes[mm].size());

            // Traverse source z-nodes
            for (ss = 0; ss < zLenSrc; ss++)
            {
                // Traverse observer z-nodes
                for (tt = 0; tt < zLenObs; tt++)
                {
                	// Set thread‑local layer state
		            mgfLocal.i = ii;
		            mgfLocal.m = mm;
		            mgfLocal.smgf.SetLayers(ii, mm);


                	// Fetch coords & row index
		            zp      = lm.z_nodes[ii][ss];
		            z       = lm.z_nodes[mm][tt];
		            idx_z   = mgfLocal.z_to_idx.at({mm, z});
		            idx_zp  = mgfLocal.z_to_idx.at({ii, zp});
		            idx_row = mgfLocal.idxpair_to_row.at({idx_z, idx_zp});

                	rho     = lm.rho_nodes[rho_idx];
					auto &row = table[idx_row];

                	// Shift everything to the right of rho_idx
					for (int j = rho_new_size - 1; j > rho_idx; --j)
						row[j] = row[j - 1];

                	// Compute the new sample
		            table_entry<N> new_entry{};
		            if (!curl) {
		                if (mgfLocal.s.sampling_method == MGF_INTEGRATE)
		                    mgfLocal.ComputeMGF_Integration(rho, z, zp, new_entry.K);
		                else
		                    mgfLocal.ComputeMGF_DCIM(rho, z, zp, new_entry.K);
		            } else {
		                if (mgfLocal.s.sampling_method == MGF_INTEGRATE)
		                    mgfLocal.ComputeCurlMGF_Integration(rho, z, zp, new_entry.K);
		                else
		                    mgfLocal.ComputeCurlMGF_DCIM (rho, z, zp, new_entry.K);
		            }

		            // Write it into the slot
		            row[rho_idx] = new_entry;

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
	int idx = 0;
	for (int ii = 0; ii < lm.z_nodes.size(); ii++)
	{
		for (int zz = 0; zz < lm.z_nodes[ii].size(); zz++)
		{
			//int idx = ii >= zz ? ii*ii + ii + zz : ii + zz*zz;
			std::pair<int,double> key = {ii, lm.z_nodes[ii][zz]};
			z_to_idx.insert(std::make_pair(key, idx));
			idx++;
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

					int idx_z = z_to_idx[{mm, z}];
					int idx_zp = z_to_idx[{ii, zp}];

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

/*! \brief Consider one layer first.*/
void MGF::AddTableMaps_z(int layer_idx, int z_idx, int z_new_size)
{

    // ====== Generate maps between z-nodes and their index in the stackup ======

    // The int pair <ii, zz> is mapped to a unique int using Szudzik's function.
    // Reference: https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way

    int map_size = z_to_idx.size();

    int idx = layer_idx >= z_new_size ? layer_idx*layer_idx + layer_idx + z_new_size : layer_idx + z_new_size*z_new_size;
    z_to_idx.insert(std::make_pair(std::pair<int,double>(layer_idx, lm.z_nodes[layer_idx][z_idx]), idx));


    // ====== Generate a map for table entries ======

    int idx_row = idxpair_to_row.size();

    // Traverse observer z-nodes
    for (int ii = 0; ii < lm.layers.size(); ii++)
    {
        for (int tt = 0; tt < lm.z_nodes[ii].size(); tt++)
        {

            double zp = lm.z_nodes[layer_idx][z_idx];
            double z = lm.z_nodes[ii][tt];

            int idx_z = z_to_idx[{ii, z}];
            int idx_zp = z_to_idx[{layer_idx, zp}];

            idxpair_to_row.insert(std::make_pair(std::make_pair(idx_z, idx_zp), idx_row));

            idx_row++;
        }
    }


    // Traverse source z-nodes
    for (int ii = 0; ii < lm.layers.size(); ii++)
    {
        for (int tt = 0; tt < lm.z_nodes[ii].size(); tt++)
        {

            double zp = lm.z_nodes[ii][tt];

            if (zp == lm.z_nodes[layer_idx][z_idx])
                continue;

            double z = lm.z_nodes[layer_idx][z_idx];

            int idx_z = z_to_idx[{layer_idx,z}];
            int idx_zp = z_to_idx[{ii,zp}];

            idxpair_to_row.insert(std::make_pair(std::make_pair(idx_z, idx_zp), idx_row));

            idx_row++;
        }
    }



    int map_size_new = z_to_idx.size();

    if (map_size_new != map_size)
        std::cout << "MGF::AddTableMaps(): Table indexing have been updated due to adaptive interpolation." << std::endl;


    return;

}

/*! \brief Consider one layer first.*/
void MGF::AddTableMaps_rho()
{

    // ====== Update the map for rho-nodes ======

    rho_to_idx.clear();

    for (int ii = 0; ii < lm.rho_nodes.size(); ii++)
        rho_to_idx.insert(std::make_pair(lm.rho_nodes[ii], ii));


    return;

}





/*! \brief Function to retrieve the nearest tabulated row for a given z-zp pair.*/
int MGF::GetRow(double z, double zp)
{

	if (!initialized)
	{
		throw std::logic_error("[ERROR] MGF:GetRow(): Table has not been generated. Call MGF::Initialize() first.");
	}


	// ====== Locate the index of the z-node nearest to z using observe layer ======
	auto z_below = z_to_idx.lower_bound({m, z});
	auto z_above = z_to_idx.upper_bound({m, z});

	int idx_z;
	double z_found;

	if (std::abs(z - z_below->first.second) < 1.0e-15)
	{
		idx_z = z_below->second;
		z_found = z_below->first.second;
	}
	else
	{
		z_below--;

		if (std::abs(z - z_below->first.second) < std::abs(z - z_above->first.second))
		{
			idx_z = z_below->second;
			z_found = z_below->first.second;
		}
		else
		{
			idx_z = z_above->second;
			z_found = z_above->first.second;
		}
	}


	// ====== Locate the index of the z-node nearest to zp using source layer  ======

	auto zp_below = z_to_idx.lower_bound({i, zp});
	auto zp_above = z_to_idx.upper_bound({i, zp});

	int idx_zp;
	double zp_found;

	if (std::abs(zp - zp_below->first.second) < 1.0e-15)
	{
		idx_zp = zp_below->second;
		zp_found = zp_below->first.second;
	}
	else
	{
		zp_below--;

		if (std::abs(zp - zp_below->first.second) < std::abs(zp - zp_above->first.second))
		{
			idx_zp = zp_below->second;
			zp_found = zp_below->first.second;
		}
		else
		{
			idx_zp = zp_above->second;
			zp_found = zp_above->first.second;
		}
	}
	

	// ====== Retrive the images at the above index pair ======

	int idx_row = idxpair_to_row[std::make_pair(idx_z, idx_zp)];

	return idx_row;

}

/*! \brief Function to retrieve the interpolation stencil for a given z and zp.*/
void MGF::GetStencil_z(double z, std::vector<int> &z_idx_stencil, std::vector<double> &z_stencil)
{

    if (!initialized)
    {
        throw std::logic_error("[ERROR] MGF:GetRow(): Table has not been generated. Call MGF::Initialize() first.");
    }
    

    // ====== Locate the index of the z-node nearest to z ======

    int layer_idx = lm.FindLayer(z);

    std::vector<double> nodes = lm.z_nodes[layer_idx];

    if (nodes.size() < s.order_z + 1)
        z_stencil = nodes;
    else
    {
        auto it = std::lower_bound(nodes.begin(), nodes.end(), z);

        int left = (it - nodes.begin()) - 1;  // iterator arithmetic
        int right = it - nodes.begin();

        for (int i = 0; i < s.order_z + 1; ++i)
        {
            if (left < 0)
            {
                z_stencil.push_back(nodes[right++]); // If 'left' is out of bounds, select from the right.
                continue;
            }
            if (right >= nodes.size())
            {
                z_stencil.push_back(nodes[left--]); // If 'right' is out of bounds, select from the left.
                continue;
            }

            // Compare the absolute difference of the values at the 'left' and 'right' pointers to the target.
            if (std::abs(z - nodes[left]) < std::abs(z - nodes[right])) {
                z_stencil.push_back(nodes[left--]); // If left is closer, select it.
            } else {
                z_stencil.push_back(nodes[right++]); // Otherwise, select the right.
            }

        }

        std::sort(z_stencil.begin(), z_stencil.end());
    }


    for (int i = 0; i < z_stencil.size(); i++)
    {
        z_idx_stencil.push_back(z_to_idx[{layer_idx, z_stencil[i]}]);
    }



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


