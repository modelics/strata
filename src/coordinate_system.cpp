/****************************** coordinate_system.cpp *******************************

 * Fixes compilation in debug, this declaration is needed
 * https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
 * Author: Damian Marek
 * Created on: Jan 10, 2019

 ***********************************************************************************/

#include "coordinate_system.hpp"

constexpr double vect3D::*vect3D::p[];


