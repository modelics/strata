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

/************************ constants.hpp ************************

 * Mathematical and other constants.
 *
 * Author: Shashwat Sharma

 ***************************************************************/


#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

namespace strata
{

	const std::complex<double> J (0.0, 1.0);
	const std::complex<double> I (0.0, 1.0);
	const double PI = M_PI;

	const double eps0 = 8.854187817e-12;	        /// Vacuum electrical permittivity
	const double mu0 = 4.0*PI*1.0e-7;			    /// Vacuum magnetic permeability
	const double c0 = 1.0/std::sqrt(eps0*mu0);		/// Vacuum speed of light
	const double eta0 = std::sqrt(mu0/eps0);		/// Vacuum wave impedance

	const double EulerGamma = 0.5772156649015328606065120900824024310421; // Euler–Mascheroni constant

	const double inf = std::numeric_limits<double>::infinity();
	
}


#endif
