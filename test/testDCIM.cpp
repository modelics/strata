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

// -----------------------------------------------------------

// The Michalski-Zheng formulation-C MGF is computed using
// the two-level DCIM for the example in Fig. 3 of
// Ling, Jin, IEEE Microw. Guided Wave Lett., 2000.
// Important note: typically, one would need surface wave pole
// extraction at high frequencies, where DCIM loses accuracy.
// Surface wave pole extraction has not yet been implemented.

// -----------------------------------------------------------

#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "MGF.hpp"


int main(int argc, char** argv)
{

	std::cout << "===========================" << std::endl;
	std::cout << "TestDCIM()" << std::endl;
	std::cout << "===========================" << std::endl;

    std::string tech_file, out_file;
	
	if (argc < 2)
	{
		throw std::runtime_error("[ERROR] Layer file was not provided.");
	}
	else if (argc < 3)
	{
		tech_file = argv[1];
		out_file = "MGFdata_DCIM.txt";
	}
	else
	{
		tech_file = argv[1];
		out_file = argv[2];
	}


	// ====== Stackup ======

	// Create a layer management object and parse the layer file
	LayerManager lm;
	lm.ProcessTechFile(tech_file, 1.0e-9);

	// Set the analysis frequency and wave number
	double f = 30.0e9;
	double omega = 2.0*M_PI*f;

	// Some useful constants are provided via the Strata namespace
	double k0 = omega*std::sqrt(strata::eps0*strata::mu0);
	double lambda0 = 2.0*M_PI/k0;

	// Precompute frequency-dependent layer data
	lm.ProcessLayers(f);

	// Print layer data to the terminal for verification
	lm.PrintLayerData(NULL, true);


	// ====== Set up the source and observation points ======

	// Set the source and observation (obs) points
	// For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 10 wavelengths away from the source point
	
	double x_src = 0.0, y_src = 0.0, z_src = 0.4e-3;
	double y_obs = 0.0, z_obs = 1.4e-3;
	
	int Nx = 500; // Number of points in the sweep
	double x_obs_min = std::abs(1.6e-4*lambda0);
	double x_obs_max = std::abs(1.6e1*lambda0);
		
	// We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
	std::vector<double> x_vec;
	strata::logspace(std::log10(x_obs_min), std::log10(x_obs_max), Nx, x_vec);


	// ====== Set z nodes for computing the complex images ======

	// Complex images are generated for each possible pair of z nodes.
	// When computing the MGF for arbitrary source and observation coordinates:
	// - along the z axis, the nearest pair of source and observation z nodes is chosen (0-order interpolation).
	
	// In this example, only two points are needed along z
	std::vector<double> z_nodes = {z_src, z_obs};

	// Provide the layer manager with the z nodes
	lm.ClearNodes_z();
	lm.InsertNodes_z(z_nodes);


	// ====== Initialize the MGF class ======
	
	// This class stores all the settings we want to use
	MGF_settings s;
	
	// This class is the "engine" which will compute the MGF
	MGF mgf;

	// Tell the class that we want to use the DCIM.
	s.method = MGF_DCIM;

	// One, two, and three levels of DCIM are supported by Strata
	// s.DCIM_method = DCIM_ONE_LEVEL;
	s.DCIM_method = DCIM_TWO_LEVEL;
	// s.DCIM_method = DCIM_THREE_LEVEL;

	// One can control the tolerances for the GPOF exponential fitting algorithm
	s.tol_svd = 1.0e-4;
	s.tol_eig = 1.0e-16;

	// One can also set an upper bound on the number of images to keep. A negative number means no bound is set.
	s.max_num_images = -1;

	// Initialize the class with the given stackup and chosen settings.
	// The interpolation table will be generated at this point, so the initialization step could take a while.
	mgf.Initialize(f, lm, s);


	// ====== Compute the MGF ======

	// Create an output file where the MGF will be exported for post-processing
	std::ofstream outfile(out_file);

	// For post-processing, we'll store the frequency and positions along the z axis in the header
	outfile << "Frequency: " << f << " Hz" << std::endl;
	outfile << "z_src: " << z_src << " m" << std::endl;
	outfile << "z_obs: " << z_obs << " m" << std::endl;

	// The lateral separation between source and observation points (rho) and all the MGF components is tabulated for each observation point
	outfile << "\nrho Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

	// The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
	// The working source and observation layers can be found with the FindLayer() method of the layer manager.
	// In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
	int i = lm.FindLayer(z_src);
	int m = lm.FindLayer(z_obs);
	mgf.SetLayers(i, m); // Source first, observation second
	
	for (int ii = 0; ii < Nx; ii++)
	{

		double x_obs = x_vec[ii];

		// In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
		double x_diff = x_obs - x_src;
		double y_diff = y_obs - y_src;

		// The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
		// In addition, there is a scalar component which is just a complex number.
		std::array<std::complex<double>, 9> G_dyadic;
		std::complex<double> G_phi;

		// Compute the MGF
		mgf.ComputeMGF(x_diff, y_diff, z_obs, z_src, G_dyadic, G_phi);

		// The dyadic MGF components are now stored in G_dyadic, while the scalar MGF is stored in G_phi.


		// ====== Optional modifications to the output ======

		// With reference to Michalski, Zheng, TAP 1990, 38 (3), equations (50)--(53), the MGF we computed above ** does include ** the cos and sin pre-factors. However, in literature, the MGF is often reported without these pre-factors. Therefore, for the purpose of this example, and to make direct comparisons to data from literature, those prefactors are cancelled out below. This section of code would not be needed in an actual MoM-type application.

		std::complex<double> zeta = std::atan2(y_diff, x_diff);
		std::complex<double> cos_term = std::cos(zeta);
		std::complex<double> sin_term = std::sin(zeta);

		if (std::abs(cos_term) > 0.0)
		{
			G_dyadic[2] /= cos_term;
			G_dyadic[6] /= cos_term;
		}
		if (std::abs(sin_term) > 0.0)
		{
			G_dyadic[5] /= sin_term;
			G_dyadic[7] /= sin_term;
		}

		
		// ====== Export data to the output text file ======

		// Lateral separation
		double rho = std::sqrt( std::pow(x_diff, 2) + std::pow(y_diff, 2) );
		outfile << rho << " " << 
				std::abs(G_dyadic[0]) << " " << std::abs(G_dyadic[1]) << " " << std::abs(G_dyadic[2]) << " " << 
				std::abs(G_dyadic[3]) << " " << std::abs(G_dyadic[4]) << " " << std::abs(G_dyadic[5]) << " " << 
				std::abs(G_dyadic[6]) << " " << std::abs(G_dyadic[7]) << " " << std::abs(G_dyadic[8]) << " " << std::abs(G_phi) << std::endl;
		
	}
		
	return 0;

}






