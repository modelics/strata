/****************************** coordinate_system.hpp *****************************

 * Classes for convenient representation of 3D coordinate systems.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 20, 2018

 **********************************************************************************/


#ifndef COORD_H
#define COORD_H


#include <iostream>
#include <cmath>
#include <sstream>


// Useful reference: http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html


/*! \brief Class to store and manipulate points or position vectors in 3D, in Cartesian coordinates.*/
class vect3D
{
public:

	double x = 0.0, y = 0.0, z = 0.0;
	

	// Constructors
	vect3D() {};
	vect3D(double x, double y, double z): x{x}, y{y}, z{z} {}

	
    // Assignment
    // vect3D& operator=(const vect3D &RHS) { x = RHS.x; y = RHS.y; z = RHS.z; return *this; };

	
    // Addition and subtraction
    vect3D& operator+=(const vect3D &RHS) { x += RHS.x; y += RHS.y; z += RHS.z; return *this; };
    vect3D& operator-=(const vect3D &RHS) { x -= RHS.x; y -= RHS.y; z -= RHS.z; return *this; };
    
    vect3D& operator+=(const double &RHS) { x += RHS; y += RHS; z += RHS; return *this; };
    vect3D& operator-=(const double &RHS) { x -= RHS; y -= RHS; z -= RHS; return *this; };

    vect3D operator+(const vect3D &RHS) { return vect3D(*this) += RHS; };
    vect3D operator-(const vect3D &RHS) { return vect3D(*this) -= RHS; };

    vect3D operator+(const double &RHS) { return vect3D(*this) += RHS; };
    vect3D operator-(const double &RHS) { return vect3D(*this) -= RHS; };

	vect3D operator+() const { return vect3D(x, y, z); };
    vect3D operator-() const { return vect3D(-x, -y, -z); };

	// friend vect3D operator+(vect3D const &LHS, vect3D const &RHS) { return LHS+RHS; };
    // friend vect3D operator-(vect3D const &LHS, vect3D const &RHS) { return LHS-RHS; };

	
    // Scalar product and division
    vect3D& operator*=(const double &RHS) { x *= RHS; y *= RHS; z *= RHS; return *this; };
    vect3D& operator/=(const double &RHS) { x /= RHS; y /= RHS; z /= RHS; return *this; };

    vect3D operator*(const double &RHS) { return vect3D(*this) *= RHS; };
    vect3D operator/(const double &RHS) { return vect3D(*this) /= RHS; };

    friend vect3D operator*(double c, vect3D &vec) { return vec*c; };
    friend vect3D operator/(double c, vect3D &vec) { return vec/c; };

	
    // Comparisons
    bool operator==(const vect3D &RHS) { return ((std::abs(x - RHS.x) <= std::abs(x)*tol) && (std::abs(y - RHS.y) <= std::abs(y)*tol) && (std::abs(z - RHS.z) <= std::abs(z)*tol)); };
    bool operator!=(const vect3D &RHS) { return !(*this == RHS); };
    
    bool operator>=(const vect3D &RHS) { return ((x >= RHS.x) && (y >= RHS.y) && (z >= RHS.z)); };
    bool operator<=(const vect3D &RHS) { return ((x <= RHS.x) && (y <= RHS.y) && (z <= RHS.z)); };
    bool operator>(const vect3D &RHS) { return !(*this <= RHS); };
    bool operator<(const vect3D &RHS) { return !(*this >= RHS); };


    // Euclidean norm
	// double norm2() const { return std::hypot(std::hypot(x, y), z); };
	double norm2() const { return std::sqrt(x*x + y*y + z*z); };
	friend double norm2(vect3D const &a) { return a.norm2(); };

	
    // Dot product
    friend double dot(vect3D const &a, vect3D const &b) { return a.x*b.x + a.y*b.y + a.z*b.z; };

    // Cross product
    friend vect3D cross(vect3D const &a, vect3D const &b) { return {(a.y*b.z - a.z*b.y), (a.z*b.x - a.x*b.z), (a.x*b.y - a.y*b.x)}; };

	
    // Print to stream
    friend std::ostream& operator<<(std::ostream &stream, vect3D const &p) { return stream << std::scientific << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::flush; };

	// Convert to string
    friend std::string to_string(vect3D const &p) { std::ostringstream out; out << p; return out.str(); };

	
	// Function to explicitly return coordinates as an array of doubles, if needed
	void DoubleArray(double *v) { v[0] = x; v[1] = y; v[2] = z; return; };

	
    // To access coordinates using square bracket notation, for convenience
	double operator[](int i) const { return this->*p[i]; };
    double& operator[](int i) { return this->*p[i]; };

	
private:
    static constexpr double vect3D::*p[] = { &vect3D::x, &vect3D::y, &vect3D::z };
	static constexpr double tol = 1.0e-15;
	
};


#endif

