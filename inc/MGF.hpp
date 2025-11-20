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

/********************************** MGF.hpp *******************************

 * Routines for computing the multilayer Green's function (MGF) based on
 * Michalski & Zheng's Formulation-C. Currently supports direct numerical
 * integration, multilevel DCIM, and quasistatic analysis and extraction.
 *
 * Author: Shashwat Sharma
 * Created on: Apr 02, 2020

 *************************************************************************/


#ifndef MGF_H
#define MGF_H


#include <complex>
#include <vector>
#include <array>
#include <string>

#include "utility.hpp"
#include "layers.hpp"
#include "spectral_MGF.hpp"
#include "quasistatic_MGF.hpp"
#include "DCIM.hpp"
#include "constants.hpp"


const int MGF_DCIM = 1;
const int MGF_INTEGRATE = 2;
const int MGF_QUASISTATIC = 3;
const int MGF_INTERPOLATE = 4;

const int MGF_GEM = 5;
const int MGF_GHJ = 6;


/*! \brief Data structure for the tabulation of the MGF. Each element of the K vector stores an individual MGF component at a particular location in the table.*/
template<std::size_t N>
struct table_entry
{
	std::array<std::complex<double>, N> K {};
};


/*! \brief Data structure to store all MGF-related settings.*/
struct MGF_settings
{

	// ====== Basic settings ======

	int method = MGF_INTERPOLATE;
	bool extract_quasistatic = false;
	bool extract_singularities = false;
	bool extract_homogeneous = false;
	bool verbose = true;
	bool compute_curl = false;
	int curl_type = MGF_GEM;
	
	// ------ DCIM ------

	int DCIM_method = DCIM_TWO_LEVEL;

	// ------ Interpolation ------

	int order = 3;
	std::string filename, filename_curl;
	bool load_table = false;
	bool export_table = false;
	int sampling_method = MGF_INTEGRATE;

	
	// ====== Advanced settings ======

	std::vector<bool> components = {1, 1, 1, 1, 1};
	std::vector<bool> components_curl = {1, 1, 1, 1};

	// ------ Integration ------

	double tol_qse = 1.0e-4;	
	double switching_point = -1.0;
	
	// ------ DCIM ------

	double tol_svd = 1.0e-4;
	double tol_eig = 1.0e-16;
	int max_num_images = -1;

};


/*! \brief Class to manage and generate the MGF.*/
class MGF
{
public:

	// ============ Interface ============

    void Initialize(double _f, LayerManager &_lm, MGF_settings &_s);
    void SetLayers(int _i, int _m);
	void SetSingularityExtraction(bool extract_singularities);
	void SetComponents(std::vector<bool> components);
	void SetCurlComponents(std::vector<bool> components_curl);
	
	void ComputeMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G, std::complex<double> &G_phi);
	void ComputeCurlMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G);
	
	void ComputeQMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G, std::complex<double> &G_phi);
	void ComputeCurlQMGF(double x_diff, double y_diff, double z, double zp, std::array<std::complex<double>, 9> &G);
	
	std::complex<double> GetSingularityFactor(int component);

	
	// ============ Computational drivers ============

	template<std::size_t N>
	void ComputeMGF_Integration(double rho, double z, double zp, std::array<std::complex<double>, N> &G);

	template<std::size_t N>
	void ComputeCurlMGF_Integration(double rho, double z, double zp, std::array<std::complex<double>, N> &G);

	template<std::size_t N>
	void ComputeMGF_DCIM(double rho, double z, double zp, std::array<std::complex<double>, N> &G);

	template<std::size_t N>
	void ComputeCurlMGF_DCIM(double rho, double z, double zp, std::array<std::complex<double>, N> &G);

	template<std::size_t N>
	void ComputeMGF_Interpolation(double rho, double z, double zp, std::array<std::complex<double>, N> &G, std::vector<std::vector<table_entry<N>>> &table, std::vector<bool> &components);

	void ComputeSingularityFactors(double x_diff, double y_diff, double z, double zp);
	void ComputeSingularityFactors();
	void ComputeHomogeneousFactors();


	// ============ Computational helpers ============

	bool UseQuasistaticOnly(double rho, double z, double zp);
	
	template<std::size_t N>
	void TabulateMGF(std::vector<std::vector<table_entry<N>>> &table, bool curl = false);
	void GenerateTableMaps();
	int GetRow(double z, double zp);
	std::vector<int> GetColumns(double rho);

	template<std::size_t N>
	int LoadTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename);
	template<std::size_t N>
	int ExportTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename);
	
	
	// ============ Storage ============

	MGF_settings s;
	LayerManager lm;
	SpectralMGF smgf;
	QuasistaticMGF qmgf;
	DCIM dcim;

	double f, omega;
	int i, m;
	std::complex<double> cos_term, sin_term;
	std::complex<double> cos2_term, sin2_term;
	std::array<std::complex<double>, 9> F;
	std::complex<double> F_phi;

	// ------ MGF tabulation ------
	
	std::vector<std::vector<table_entry<5>>> MGF_table;
	std::vector<std::vector<table_entry<4>>> CurlMGF_table;

	// Maps
	std::map<std::pair<int, double>, int> z_to_idx;
	std::map<double, int> rho_to_idx;
	std::map<std::pair<int, int>, int> idxpair_to_row;
	
	// ------ Switchboard ------
	
	bool initialized = false;
	bool layers_set = false;
	bool singularity_factors_computed = false;
		

};


#endif


