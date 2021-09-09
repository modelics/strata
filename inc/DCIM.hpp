/************************** DCIM.hpp ***************************

 * Functionality to compute the spatial domain MGF using the 
 * one-, two-, or three-level DCIM with GPOF.
 * 
 * References:
 *
 * Two-level DCIM:
 * Aksun, "A Robust Approach for the Derivation of Closed-Form
 * Green's Functions", IEEE-MTT, 1996
 * Aksun, Dural, "Clarification of Issues on the Closed-Form
 * Green's Functions in Stratified Media", IEEE-TAP, 2005
 *
 * Three-level DCIM:
 * Alparslan, Aksun, Michalski, "Closed-Form Green’s Functions
 * in Planar Layered Media for All Ranges and Materials",
 * IEEE-MTT, 2010
 *
 * Modifications to the DCIM path:
 * Shuley, Boix, Medina, Horno, "On the Fast Approximation of
 * Green’s Functions in MPIE Formulations for Planar Layered
 * Media", IEEE-MTT, 2002
 * Kipp, Chan, "Complex Image Method for Sources in Bounded
 * Regions of Multilayer Structures", IEEE-MTT, 1994
 *
 * Author: Shashwat Sharma
 * Date Created: Mar 31, 2020

 ***************************************************************/


#ifndef DCIM_H
#define DCIM_H


#include <complex>
#include <vector>
#include <map>

#include "layers.hpp"
#include "spectral_MGF.hpp"


const int DCIM_ONE_LEVEL = 1;
const int DCIM_TWO_LEVEL = 2;
const int DCIM_THREE_LEVEL = 3;


/*! \brief Data structure to store DCIM settings.*/
struct DCIM_settings
{
	std::vector<bool> components = {1, 1, 1, 1, 1};
	std::vector<bool> components_curl = {1, 1, 1, 1};
    int method = DCIM_TWO_LEVEL;
	bool extract_quasistatic = false;
	bool extract_homogeneous = false;

	bool compute_curl = false;
	bool curl_GEM = false;

	// GPOF settings
	double tol_svd = 1.0e-4;
	double tol_eig = 1.0e-16;
	int max_num_images = -1;
};


/*! \brief Data structure to store computed DCIM images.*/
struct DCIM_images
{
	std::complex<double> k;
	std::vector<std::complex<double>> a, alpha;
};


/*! \brief Data structure to store lateral and axial wave number samples for DCIM.*/
struct DCIM_sample_points
{
	std::complex<double> k;
	double T01, T02, T03;
	std::vector<double> t1, t2, t3;
	std::vector<std::complex<double>> krho_path1, krho_path2, krho_path3;
	std::vector<std::complex<double>> kz_path1, kz_path2, kz_path3;
};


/*! \brief Class to generate and retrieve DCIM images.*/
class DCIM
{
 public:

	// ====== Interface ======

	void GenerateImages(LayerManager *_lm, double _f, DCIM_settings _s);
	std::vector<DCIM_images>& GetImages(double z, double zp, bool curl_mode = false);
	

	// ====== Two-level DCIM ======

	void GenerateSamplePoints_TwoLevel(DCIM_sample_points &sp, int N1 = 100, int N2 = 100);
	void GenerateImages_TwoLevel(double z, double zp, std::vector<DCIM_images> &im, bool curl_mode = false);


	// ====== One-level DCIM ======

	void GenerateSamplePoints_OneLevel(DCIM_sample_points &sp, int N1 = 100);
	void GenerateImages_OneLevel(double z, double zp, std::vector<DCIM_images> &im, bool curl_mode = false);


	// ====== Three-level DCIM ======

	void GenerateSamplePoints_ThreeLevel(DCIM_sample_points &sp, int N1 = 100, int N2 = 100, int N3 = 100);
	void GenerateImages_ThreeLevel(double z, double zp, std::vector<DCIM_images> &im, bool curl_mode = false);


	// ====== Computational helpers ======

	void GenerateSpectralSamples(SpectralMGF &smgf, std::vector<std::vector<std::complex<double>>> &K, std::vector<std::complex<double>> &krho, std::vector<std::complex<double>> &kz, int N, int N_components, bool curl_mode = false);
	void RunGPOF(std::vector<std::complex<double>> &y, double dt, std::vector<std::complex<double>> &a, std::vector<std::complex<double>> &alpha, double tol_svd = 1.0e-4, double tol_eig = 1.0e-16, int max_num_images = -1);


	// ====== Storage ======

	LayerManager *lm;
    double f, omega;
    DCIM_settings s;

	// Images
	bool images_generated = false;
	std::vector<std::vector<DCIM_images>> im, imc;

	// Maps
	std::map<double, int> z_to_idx;
	std::map<std::pair<int, int>, int> idxpair_to_image;
	

};


#endif
