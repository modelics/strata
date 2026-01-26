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
    bool interpolate_z = false;
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
    int order_z = 3;
    double N_lambda = 10.0;
    double adaptive_threshold = 0.01;
    bool update_z_nodes = false;
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
    void AdaptiveInterpolation();
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

    template<std::size_t N>
    void ComputeMGF_Interpolation_withZ(double rho, double z, double zp, std::array<std::complex<double>, N> &G, std::vector<std::vector<table_entry<N>>> &table, std::vector<bool> &components);

	void ComputeSingularityFactors(double x_diff, double y_diff, double z, double zp);
	void ComputeSingularityFactors();
	void ComputeHomogeneousFactors();
	
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
		
private:

	// ============ Testing & Debugging ============

    void PlotRhoNodes(std::vector<double> z_gridpoints);
    void ProcessRhoTests(double rho_spacing, double z_test, double z_src, int level, std::vector<double>& rho_tests);
    void ProcessZTests(int layer_idx, double rho_test, double z_spacing, double z_test, double z_src, int level, std::vector<double>& z_tests);
    void TestAddRhoTableRecursive(double rho_test_l1, double rho_spacing, double z_test, double z_src, int level, std::vector<double> &rho_test_l2);
    void TestAddZTableRecursive(int layer_idx, double rho_test, double z_spacing, double z_test, double z_src, int level, std::vector<double> &z_tests_l2);

	// ============ Computational helpers ============

	bool UseQuasistaticOnly(double rho, double z, double zp);
	
	template<std::size_t N>
	void TabulateMGF(std::vector<std::vector<table_entry<N>>> &table, bool curl = false);
    template<std::size_t N>
    void AppendMGFTableZ(std::vector<std::vector<table_entry<N>>> &table, int layer_idx, int z_idx, int z_new_idx, bool curl = false);
    template<std::size_t N>
    void AppendMGFTableRho(std::vector<std::vector<table_entry<N>>> &table, int rho_idx, int rho_new_idx, bool curl = false);
    void GenerateTableMaps();
    void AddTableMapsZ(int layer_idx, int z_idx, int z_new_idx);
    void AddTableMapsRho();
	int GetRow(double z, double zp);
    void GetStencilZ(double z, std::vector<int> &z_idx_stencil, std::vector<double> &z_stencil);
	std::vector<int> GetColumns(double rho);
    void UpdateZNodes();
    void UpdateRhoNodes();
    bool IsMidpointCorrect(double rho, double z_src, double z_test,  double adaptive_threshold, bool test_rho);
    void AddRhoTable(double rho_test);
    void AddZTable(double z_test, int layer);

	template<std::size_t N>
	int LoadTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename);
	template<std::size_t N>
	int ExportTable(std::vector<std::vector<table_entry<N>>> &table, std::string filename);

};

#endif


