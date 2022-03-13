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

// -----------------------------------------------------------

// The computation of the multilayer Green's function (MGF)
// in spectral domain is demonstrated in this example.
// The spectral MGF is computed for the example in Fig. 3
// in Yuan, Sarkar, Salazar-Palma, IEEE Trans. Microw. Theory 
// Tech., vol. 54, no. 3, 2006.

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
	std::cout << "TestSpectralMGF()" << std::endl;
	std::cout << "===========================" << std::endl;

    std::string tech_file, out_file;
	
	if (argc < 2)
	{
		throw std::runtime_error("[ERROR] Layer file was not provided.");
	}
	else if (argc < 3)
	{
		tech_file = argv[1];
		out_file = "MGFdata.txt";
	}
	else
	{
		tech_file = argv[1];
		out_file = argv[2];
	}


	// ====== Stackup ======

	// Create a layer management object and parse the layer file
	LayerManager lm;
	lm.ProcessTechFile(tech_file);

	// Set the analysis frequency and wave number
	double f = 30.0e9;
	double omega = 2.0*M_PI*f;

	// Some useful constants are provided via the Strata namespace
	double eps0 = strata::eps0;
	double mu0 = strata::mu0;
	std::complex<double> J = strata::J; // imaginary unit

	double k0 = omega*std::sqrt(eps0*mu0);
	double lambda0 = 2.0*M_PI/k0;

	// Precompute frequency-dependent layer data
	lm.ProcessLayers(f);

	// Print layer data to the terminal for verification
	lm.PrintLayerData(NULL, true);


	// ====== Set up the source and observation points ======

	// Set the source and observation (obs) points
	// For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 10 wavelengths away from the source point
	
	double x_src = 0.0, y_src = 0.0, z_src = 1.0e-3;
	double y_obs = 0.0, z_obs = 1.0e-3;

	// The source and observation points are at the layer interface and could be seen as belonging to either layer. Let's just assign them to be in layer 0 (not necessary - Strata can handle this internally - but just for simplicity)
	int i = 0, m = 0;

	int Nt = 500; // Number of points to sample the spectral MGF along the deformed path in the complex plane
	double t_min = 0.0;
	double T0 = 15.0;
		
	// We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
	std::vector<double> t_vec;
	strata::linspace(t_min, T0, Nt, t_vec);


	// ====== Initialize the spectral MGF class ======

	// Initialize and set up the spectral MGF object
	SpectralMGF smgf;
	smgf.Initialize(&lm, f);
	smgf.SetLayers(i, m);

	smgf.SetSourcePoint(z_src);
	smgf.SetObservationPoint(z_obs);


	// ====== Compute the spectral MGF ======
	
	// Wave number in the source layer
	std::complex<double> k = lm.k[i];
	// k = k0;
	std::vector<std::complex<double>> G (Nt);

	for (int ii = 0; ii < Nt; ii++)
	{

		// The spectral MGF is computed for a given krho. We need to map t defined above to krho (see equation (10) of the Yuan paper)
		std::complex<double> kz = k0*(-J*t_vec[ii] + (1.0 - t_vec[ii]/T0));
		std::complex<double> krho = std::sqrt(std::pow(k0, 2) - std::pow(kz, 2));

		// Ensure that the samples in krho space lie in the top-right complex quadrant
		krho.real(std::abs(std::real(krho)));
		krho.imag(std::abs(std::imag(krho)));

		// Apply a scaling factor to allow comparing results to Yuan 2006, based on eqn. (4)
		std::complex<double> scaling_factor = 2.0*J*kz;///std::exp(-J*kz*std::abs(z_src + z_obs));
		// scaling_factor *= 4.0*M_PI/mu0/J/k;
		// scaling_factor = 1.0;

		smgf.SetRadialWaveNumber(krho);
		smgf.ComputeSpectralMGF();

		// We only care about the scalar potential part for this example
		// For some reason, there is an offset in the real part between the results here and those in Yuan 2006, possibly due to the static term extraction they mention
		G[ii] = smgf.GetResult(4)*scaling_factor - 1.0;

	}


	// ====== Export to file ======
	
	// Create an output file where the spectral MGF will be exported for post-processing
	std::ofstream outfile(out_file);

	// For post-processing, we'll store the frequency and positions along the z axis in the header
	outfile << "Frequency: " << f << " Hz" << std::endl;
	outfile << "z_src: " << z_src << " m" << std::endl;
	outfile << "z_obs: " << z_obs << " m" << std::endl;

	outfile << "\nt Re(Gphi) Im(Gphi)" << std::endl;

	for (int ii = 0; ii < Nt; ii++)
		outfile << t_vec[ii] << " " << std::real(G[ii]) << " " << std::imag(G[ii]) << std::endl;
		
	return 0;

}






