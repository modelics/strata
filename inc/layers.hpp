/*************************** layers.hpp ******************************
 
 * Classes to store and manage stratified background media.
 *
 * Author: Shashwat Sharma
 * Date Created: July 17, 2019

 *********************************************************************/


#ifndef LAYERS_H
#define LAYERS_H


#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


/*! \brief Class to store and manage frequency-dependent information for each layer.*/
class Layer
{
public:

	// ====== Storage ======

	int layerID;
	double zmax, zmin, h, epsr, sigma, mur, sigmamu;
	std::vector<int> objectIDs;   /// List of indices of objects that belong to this layer
	std::vector<double> obj_zmin, obj_zmax;  /// Min and max z extent of each object


	// ====== Initialization ======
	
	Layer() {};
    
    Layer(double _zmin, double _zmax, double _epsr = 1.0, double _mur = 1.0, double _sigma = 0.0, double _sigmamu = 0.0)
    {
		epsr = _epsr;
		sigma = _sigma;
		mur = _mur;
		sigmamu = _sigmamu;

		zmin = _zmin;
		zmax = _zmax;
		h = zmax - zmin;

		if (h <= 0.0)
			std::cout << "[ERROR] Layer(): A layer with invalid height (" << h << ") was added." << std::endl;

		return;
    };

};


/*! \brief Class to store and manage all layered medium information.*/
class LayerManager
{
public:

  	// ====== Interface - layer management ======

	void ProcessTechFile(std::string tech_file, double units = 1.0e-9);
	void ProcessTechFile_yaml(std::string tech_file, double units = 1.0e-9);
	void ProcessTechFile_tech(std::string tech_file, double units = 1.0e-9);
	
	void AddLayer(double zmin, double zmax, double epsr, double mur, double sigma, double sigmamu);
	void SetHalfspaces(double _epsr_top, double _mur_top, double _sigma_top, double _epsr_bot, double _mur_bot, double _sigma_bot, bool _isPEC_top, bool _isPEC_bot);
	
	void ProcessLayers(double f);

	int FindLayer(double z);

	void PrintLayerData(std::ofstream *out_file = NULL, bool print_to_terminal = true);
	
  
	// ====== Interface - node management ======

	void InsertNodes_z(std::vector<double> &nodes);
	void InsertNodes_z(std::vector<double> &nodes, int idx_layer);
	void PruneNodes_z(double dist);
	void PruneNodes_z(double dist, int idx_layer);
	void ClearNodes_z();
	
	void InsertNodes_rho(std::vector<double> &nodes);
	void ClearNodes_rho();

	void PrintNodeData_z(std::ofstream *out_file, bool print_to_terminal);
	void PrintNodeData_rho(std::ofstream *out_file, bool print_to_terminal);
	

	// ====== Helpers - layer management ======

	void MergeLayersWithSameMaterial(double tol = 1.0e-15);
	double GetZmin(int idx);
	double GetZmax(int idx);
	double GetHeight(int idx);
	
  
	// ====== Storage ======

	bool layers_processed = false;
	std::vector<Layer> layers;
	std::vector<std::complex<double>> k;    /// Complex wave number of each layer
	std::vector<std::complex<double>> eps;  /// Complex permittivity of each layer
	std::vector<std::complex<double>> mu;   /// Complex permeability of each layer

	std::complex<double> k_min, k_max;      /// Smallest and largest wave numbers in the structure, including bounding half spaces
	int idx_k_min, idx_k_max;               /// Layer indices of the layers that have the smallest and largest wave numbers
	
	// ------ Top and bottom halfspace properties ------
	
	double epsr_top = 1.0, mur_top = 1.0, sigma_top = 0.0;
	double epsr_bot = 1.0, mur_bot = 1.0, sigma_bot = 0.0;

	bool isPEC_top = false;
	bool isPEC_bot = false;

	std::complex<double> eps_top, mu_top, k_top;
	std::complex<double> eps_bot, mu_bot, k_bot;

	// ------ MGF (pre)computation nodes and settings ------

	// Container to store lists of nodes along z (inner vector) per layer (outer vector) on which MGF (pre)computations are to take place.
	std::vector<std::vector<double>> z_nodes;

	// Vector to store lists of nodes along rho on which MGF (pre)computations are to take place.
	std::vector<double> rho_nodes;
	
};


#endif
