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

/*************************** layers.cpp ******************************
 
 * Classes to store and manage stratified background media.
 *
 * Author: Shashwat Sharma
 * Date Created: July 17, 2019

 *********************************************************************/


#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "yaml-cpp/yaml.h"

#include "layers.hpp"
#include "constants.hpp"
#include "utility.hpp"

using namespace strata;


// ==================================================================================
// Interface - layer management
// ==================================================================================

/*! \brief Function to parse a given technology file and populate layer data.*/
void LayerManager::ProcessTechFile(std::string tech_file, double units)
{

	std::string extension;
	std::string::size_type idx = tech_file.rfind('.');

	if (idx != std::string::npos)
	    extension = tech_file.substr(idx+1);
	else
	    throw std::invalid_argument("[ERROR] LayerManager::ProcessTechFile(): Technology file does not seem to have an extension. Extension must be either .yaml or .tech.");

	if (extension == "yaml" || extension == "yml")
		ProcessTechFile_yaml(tech_file, units);
	else if (extension == "tech")
		ProcessTechFile_tech(tech_file, units);
	else
	    throw std::invalid_argument("[ERROR] LayerManager::ProcessTechFile(): Technology file has an invalid extension. Extension must be either .yaml or .tech.");
	
	return;

}


/*! \brief Function to parse a given .yaml file and populate layer data.*/
void LayerManager::ProcessTechFile_yaml(std::string tech_file, double units)
{

	layers_processed = false;
	
	YAML::Node node = YAML::LoadFile(tech_file);

	std::string unit_str = node["unit"].as<std::string> ("nm");
	if (unit_str == "m")
		units = 1.0;
	else if (unit_str == "mm")
		units = 1.0e-3;
	else if (unit_str == "um")
		units = 1.0e-6;
	else if (unit_str == "nm")
		units = 1.0e-9;
	else
		throw std::invalid_argument("[ERROR] LayerManager::ProcessTechFile_yaml(): Invalid units provided in the technology file. Units must be one of {m, mm, um, nm}");

	// Define the default top and bottom half-spaces as air
	bool _isPEC_top = false;
	bool _isPEC_bot = false;

	double _epsr_top = 1.0;
	double _mur_top = 1.0;
	
	double _epsr_bot = 1.0;
	double _mur_bot = 1.0;
	
	double _sigma_top = 0.0;
	double _sigma_bot = 0.0;


	// ====== Dielectric layers ======

	if (const YAML::Node &diel = node["dielectric_layers"])
	{
		for (YAML::const_iterator it = diel.begin(); it != diel.end(); ++it)
		{

			std::string label = it->first.as<std::string> ();
			const YAML::Node &layer = it->second;

			double zmin = layer["zmin"].as<double> ();
			double h = layer["h"].as<double> ();
			
			zmin *= units;
			h *= units;
			double zmax = zmin + h;

			double epsr = layer["epsr"].as<double> ();
			double mur = layer["mur"].as<double> (1.0);
			double sigma = layer["sigma"].as<double> (0.0);
			double sigmamu = layer["sigmamu"].as<double> (0.0);

			// Add this layer to the full layer-set
			AddLayer(zmin, zmax, epsr, mur, sigma, sigmamu);
			
		}
	}
	else
	{
		throw std::runtime_error("[ERROR] LayerManager::ProcessTechFile_yaml(): No dielectric layers were provided in the technology file.");
	}


	// ====== Top halfspace ======

	if (const YAML::Node &top = node["top_halfspace"])
	{
		_epsr_top = top["epsr"].as<double> (1.0);
		_mur_top = top["mur"].as<double> (1.0);
		_sigma_top = top["sigma"].as<double> (0.0);

		if (_sigma_top < 0)
		{
			_isPEC_top = true;
			_sigma_top = 0.0;
		}
	}


	// ====== Bottom halfspace ======

	if (const YAML::Node &bot = node["bottom_halfspace"])
	{
		_epsr_bot = bot["epsr"].as<double> (1.0);
		_mur_bot = bot["mur"].as<double> (1.0);
		_sigma_bot = bot["sigma"].as<double> (0.0);

		if (_sigma_bot < 0)
		{
			_isPEC_bot = true;
			_sigma_bot = 0.0;
		}
	}

	SetHalfspaces(_epsr_top, _mur_top, _sigma_top, _epsr_bot, _mur_bot, _sigma_bot, _isPEC_top, _isPEC_bot);
	
	
	return;

}


/*! \brief Function to parse a given .tech file and populate layer data.*/
void LayerManager::ProcessTechFile_tech(std::string tech_file, double units)
{

	layers_processed = false;
	
	std::ifstream tech_stream;
	tech_stream.open(tech_file);

	if (!tech_stream.good())
		std::cout << "[ERROR] LayerManager::ProcessTechFile(): Tech file could not be found or opened." << std::endl;

	std::string label;
	double zmin, zmax, h;
	double epsr, mur, sigma, sigmamu = 0.0;

	// Define the default top and bottom half-spaces as air
	bool _isPEC_top = false;
	bool _isPEC_bot = false;

	double _epsr_top = 1.0;
	double _mur_top = 1.0;
	
	double _epsr_bot = 1.0;
	double _mur_bot = 1.0;
	
	double _sigma_top = 0.0;
	double _sigma_bot = 0.0;

	std::string layer_type;
        
	// Parse
	while (!tech_stream.eof())
	{
		
		std::streampos current_line = tech_stream.tellg(); // Remember where this line starts
    
		std::string line;
		getline(tech_stream, line); // Now the current position is the end of this line
		
		// Get rid of the trailing newline character
		line = line.substr(0, line.length() - 1);
                       
		if (!line[0] || (line[0] == '/' && line[1] == '/')) // Blank line or comment
			continue;
		else if (line == "Dielectric" || line == "dielectric" || line == "metal" || line == "Metal")
		{
			layer_type = line;
			continue;
		}

		tech_stream.seekg(current_line); // Go back to the start of this line using remembered position
                
		// First column: layer ID
		int temp;
		if (tech_stream >> temp){} else break;
                
		// Second column: label
		tech_stream >> label;
                
		// Third column: z_min
		tech_stream >> zmin;

		// Fourth column: h
		tech_stream >> h;

		zmin *= units;
		h *= units;
		zmax = zmin + h;
                
		// Fifth column: epsr
		tech_stream >> epsr;

		// Sixth column: mur
		tech_stream >> mur;

		// Seventh column: sigma (S/m)
		tech_stream >> sigma;


		// Check if this layer corresponds to the top or bottom half-space. In either case, update the relevant layer-set boundary and move on to the next entry.
		if (label == "top" || label == "Top")
		{
			_epsr_top = epsr;
			_mur_top = mur;
			_sigma_top = sigma;

			// A negative conductivity is used to specify a PEC half-space
			if (sigma < 0)
			{
				_isPEC_top = true;
				_sigma_top = 0.0;
			}

			continue;
		}
		else if (label == "bot" || label == "Bot")
		{
			_epsr_bot = epsr;
			_mur_bot = mur;
			_sigma_bot = sigma;

			// A negative conductivity is used to specify a PEC half-space
			if (sigma < 0)
			{
				_isPEC_bot = true;
				_sigma_bot = 0.0;
			}

			continue;
		}


		if (layer_type == "Dielectric" || layer_type == "dielectric")
		{
			// Add this layer to the full layer-set
			AddLayer(zmin, zmax, epsr, mur, sigma, sigmamu);
		}

	}

	SetHalfspaces(_epsr_top, _mur_top, _sigma_top, _epsr_bot, _mur_bot, _sigma_bot, _isPEC_top, _isPEC_bot);
  
	tech_stream.close();

	
	return;
	
}


/*! \brief Function to add a new layer to the layer set. It is the user's responsibility to make sure that the z_min and z_max do not overlap with the z-extent of any other layer in the set.*/
void LayerManager::AddLayer(double zmin, double zmax, double epsr, double mur, double sigma, double sigmamu)
{

	if (zmax - zmin <= 0.0)
	{
		std::cout << "[WARNING] LayerManager::AddLayer(): Tried to add a layer with 0 or negative height - ignoring it and moving on." << std::endl;
	    return;
	}

	layers_processed = false;
	
	// Layers are stored in order along the z axis, from top to bottom, so the new layer needs to be inserted at the correct location. This is inferred by comparing the central z coordinates of the new layer with existing layers.

	int position = layers.size();
	double znew = (zmax + zmin)/2.0;

	for (int ii = 0; ii < layers.size(); ii++)
	{
		double z = (layers[ii].zmax + layers[ii].zmin)/2.0;

		if (znew > z)
	    {
			position = ii;
			break;
	    }
	}

	// Create a new layer
	Layer layer (zmin, zmax, epsr, mur, sigma, sigmamu);
	layers.insert(layers.begin() + position, layer);
	
	return;
	
}


/*! \brief Function to set the material properties of the upper and lower halfspaces that bound the set of layers.*/
void LayerManager::SetHalfspaces(double _epsr_top, double _mur_top, double _sigma_top, double _epsr_bot, double _mur_bot, double _sigma_bot, bool _isPEC_top, bool _isPEC_bot)
{

	layers_processed = false;
	
	epsr_top = _epsr_top;
	mur_top = _mur_top;
	sigma_top = _sigma_top;
	epsr_bot = _epsr_bot;
	mur_bot = _mur_bot;
	sigma_bot = _sigma_bot;

	isPEC_top = _isPEC_top;
	isPEC_bot = _isPEC_bot;

	return;

}


/*! \brief Function to compute frequency-dependent properties of all layers. Must be called whenever a change is made to the layer set.*/
void LayerManager::ProcessLayers(double f)
{

	auto epsilon_complex = [] (double eps0, double epsr, double sigma, double omega)
		{
			if (omega > 0.0)
				return std::complex<double> ((epsr*eps0), (-sigma/omega));
			else
				return std::complex<double> ((epsr*eps0), 0.0);
		};
	
	// Make sure that any consecutive layers with the same material properties are merged
	MergeLayersWithSameMaterial();
  
	double omega = 2.0*M_PI*f;

	// ------ Half space parameters ------

	eps_top = epsilon_complex(eps0, epsr_top, sigma_top, omega);
	eps_bot = epsilon_complex(eps0, epsr_bot, sigma_bot, omega);
	
	mu_top = std::complex<double> ((mur_top * mu0), 0.0);
	mu_bot = std::complex<double> ((mur_bot * mu0), 0.0);

	k_top = omega*std::sqrt(eps_top*mu_top);
	k_bot = omega*std::sqrt(eps_bot*mu_bot);

	// ------ Layer parameters ------
	
	eps.clear();
	mu.clear();
	k.clear();

	eps.resize(layers.size());
	mu.resize(layers.size());
	k.resize(layers.size());

	k_min = 1.0e300;
	k_max = 0.0;

	for (int ii = 0; ii < layers.size(); ii++)
	{

		layers[ii].layerID = ii;

		eps[ii] = epsilon_complex(eps0, layers[ii].epsr, layers[ii].sigma, omega);
		mu[ii] = std::complex<double> ((layers[ii].mur * mu0), 0.0);

		k[ii] = omega*std::sqrt(eps[ii]*mu[ii]);

		if (std::abs(k[ii]) < std::abs(k_min))
		{
			k_min = k[ii];
			idx_k_min = ii;
		}
		
		if (std::abs(k[ii]) > std::abs(k_max))
		{
			k_max = k[ii];
			idx_k_max = ii;
		}
		
	}
	
	// if (std::abs(k_top) < std::abs(k_min))
	// {
	// 	k_min = k_top;
	// 	idx_k_min = -1;
	// }
	
	// if (std::abs(k_bot) < std::abs(k_min))
	// {
	// 	k_min = k_bot;
	// 	idx_k_min = layers.size();
	// }
	
	// if (std::abs(k_top) > std::abs(k_max))
	// {
	// 	k_max = k_top;
	// 	idx_k_max = -1;
	// }
	
	// if (std::abs(k_bot) > std::abs(k_max))
	// {
	// 	k_max = k_bot;
	// 	idx_k_max = layers.size();
	// }

	layers_processed = true;

	return;
	
}


/*! \brief Function to find the layer within which a particular z-coordinate exists. Returns -1 if no eligible layer was found.*/
int LayerManager::FindLayer(double z)
{

	double tol = (layers[0].zmax - layers.back().zmin)*1.0e-15;

	if (z > layers[0].zmax)
		return -1;
	if (z < layers.back().zmin)
		return layers.size();

	for (int ii = 0; ii < layers.size(); ii++)
	{
		if (layers[ii].zmax - z > tol && z - layers[ii].zmin > tol)
			// The point is in layer ii, and not at an interface
			return ii;
		else if (std::abs(layers[ii].zmax - z) <= tol)
			// The point is at the top interface of this layer
			return ii;
		else if (ii == (int)(layers.size() - 1))
			// The point must be at the bottom interface of the bottom layer
			return ii;
	}

	throw std::out_of_range("[ERROR] LayerManager::FindLayer(): Could not locate z-point " + std::to_string(z) + ".");
	
}


/*! \brief Function to print layer information for debugging and logging purposes.*/
void LayerManager::PrintLayerData(std::ofstream *out_file, bool print_to_terminal)
{

	std::string message;

	message = "\n=============================================\n";
	message += "Layer data\n";
	message += "=============================================\n\n";

	message += "Number of layers: " + std::to_string(layers.size()) + "\n\n";

	message += "------------------------------\n";
    message += "Top half-space:\n";
    message += "Is PEC: " + std::to_string(isPEC_top) + "\n";
	message += "Relative permittivity: " + std::to_string(epsr_top) + "\n";
	message += "Relative permeability: " + std::to_string(mur_top) + "\n";
	message += "Electrical conductivity: " + std::to_string(sigma_top) + "\n";
	if (layers_processed)
		message += "Wave number: (" + std::to_string(k_top.real()) + ") + j(" + std::to_string(k_top.imag()) + ")\n";
		
	for (int ii = 0; ii < layers.size(); ii++)
	{
		message += "------------------------------\n";
		message += "Layer ID: " + std::to_string(layers[ii].layerID) + "\n";
		message += "zmax: " + std::to_string(layers[ii].zmax) + "\n";
		message += "zmin: " + std::to_string(layers[ii].zmin) + "\n";
		message += "Height: " + std::to_string(layers[ii].h) + "\n";
		message += "Relative permittivity: " + std::to_string(layers[ii].epsr) + "\n";
		message += "Relative permeability: " + std::to_string(layers[ii].mur) + "\n";
		message += "Electrical conductivity: " + std::to_string(layers[ii].sigma) + "\n";
		if (layers_processed)
			message += "Wave number: (" + std::to_string(k[ii].real()) + ") + j(" + std::to_string(k[ii].imag()) + ")\n";
	}

    message += "------------------------------\n";
	message += "Bottom half-space:\n";
    message += "Is PEC: " + std::to_string(isPEC_bot) + "\n";
	message += "Relative permittivity: " + std::to_string(epsr_bot) + "\n";
	message += "Relative permeability: " + std::to_string(mur_bot) + "\n";
	message += "Electrical conductivity: " + std::to_string(sigma_bot) + "\n";
	if (layers_processed)
		message += "Wave number: (" + std::to_string(k_bot.real()) + ") + j(" + std::to_string(k_bot.imag()) + ")\n";
	message += "------------------------------\n";
		
	message += "\n=============================================\n";

	if (print_to_terminal)
		std::cout << message << std::flush;
	if (out_file != NULL)
		(*out_file) << message << std::flush;

	return;

}


// ==================================================================================
// Interface - node management
// ==================================================================================

/*! \brief Function to insert the set of z-coordinates over which precomputations should take place. Note that this function adds the given nodes to the existing list of nodes. Call ClearNodes_z() to reset the list.*/
void LayerManager::InsertNodes_z(std::vector<double> &nodes)
{

	// The z-coordinates need to be inserted layer-wise. Convention: a node that sits right at an interface is considered to belong to the layer above it. The only exception is the top-most layer; nodes that lie on the interface between the top-most surface and the top half-space are assumed to belong to the top-most layer.

	z_nodes.resize(layers.size());

	for (int ii = 0; ii < layers.size(); ii++)
	{

		double tol = 1.0e-6*layers[ii].h;
		
		for (int jj = 0; jj < nodes.size(); jj++)
		{
			// Top-most layer
			if (ii == 0 && nodes[jj] <= layers[ii].zmax && nodes[jj] >= layers[ii].zmin)
				z_nodes[ii].push_back(nodes[jj]);

			// All other layers
			else if (nodes[jj] < layers[ii].zmax && nodes[jj] >= layers[ii].zmin)
				z_nodes[ii].push_back(nodes[jj]);
		}

		// Sort and erase duplicates
		std::sort(z_nodes[ii].begin(), z_nodes[ii].end());
		z_nodes[ii].erase(std::unique(z_nodes[ii].begin(), z_nodes[ii].end(), [tol](double a, double b) { return std::abs(a - b) <= tol; }), z_nodes[ii].end());

	}

	return;

}


/*! \brief Function to insert to a given layer the set of z-coordinates over which precomputations should take place. Note that this function adds the given nodes to the existing list of nodes. Call ClearNodes_z() to reset the list.*/
void LayerManager::InsertNodes_z(std::vector<double> &nodes, int idx_layer)
{

	double tol = 1.0e-6*layers[idx_layer].h;
	
	z_nodes.resize(layers.size());
	z_nodes[idx_layer].insert(z_nodes[idx_layer].end(), nodes.begin(), nodes.end());

	// Sort and erase duplicates
	std::sort(z_nodes[idx_layer].begin(), z_nodes[idx_layer].end());
	z_nodes[idx_layer].erase(std::unique(z_nodes[idx_layer].begin(), z_nodes[idx_layer].end(), [tol](double a, double b) { return std::abs(a - b) <= tol; }), z_nodes[idx_layer].end());

	return;

}


/*! \brief Delete z-coordinates which are within a given physical distance of each other, for all layers.*/
void LayerManager::PruneNodes_z(double dist)
{
	for (int ii = 0; ii < layers.size(); ii++)
		PruneNodes_z(dist, ii);
	return;
}


/*! \brief Delete z-coordinates which are within a given physical distance of each other, for a given layer.*/
void LayerManager::PruneNodes_z(double dist, int idx_layer)
{

	// Sort (just in case) and erase
	std::sort(z_nodes[idx_layer].begin(), z_nodes[idx_layer].end());
	z_nodes[idx_layer].erase(std::unique(z_nodes[idx_layer].begin(), z_nodes[idx_layer].end(), [dist](double a, double b) { return std::abs(a - b) <= dist; }), z_nodes[idx_layer].end());

	return;

}


/*! \brief Clear the existing list of z-nodes.*/
void LayerManager::ClearNodes_z() { z_nodes.clear(); return; }


/*! \brief Function to insert the set of rho-coordinates over which precomputations should take place. Only relevant for interpolation-based MGF approaches. Note that this function adds the given nodes to the existing list of nodes. Call ClearNodes_rho() to reset the list.*/
void LayerManager::InsertNodes_rho(std::vector<double> &nodes)
{

	// Append to the existing set of nodes
	rho_nodes.insert(rho_nodes.end(), nodes.begin(), nodes.end());

	// Sort and erase duplicates
	std::sort(rho_nodes.begin(), rho_nodes.end());
	rho_nodes.erase(std::unique(rho_nodes.begin(), rho_nodes.end() ), rho_nodes.end());
	
	return;

}


/*! \brief Clear the existing list of rho-nodes.*/
void LayerManager::ClearNodes_rho() { rho_nodes.clear(); return; }


/*! \brief Function to print z-node information for debugging and logging purposes.*/
void LayerManager::PrintNodeData_z(std::ofstream *out_file, bool print_to_terminal)
{

	std::string message;

	message = "\n=============================================\n";
	message += "Node data (z-nodes)\n";
	message += "=============================================\n\n";

	int N_nodes = 0;
	for (int ii = 0; ii < z_nodes.size(); ii++)
		N_nodes += z_nodes[ii].size();
	
	message += "Number of nodes: " + std::to_string(N_nodes) + "\n";

	message += "\n-----------------------------------------\n";
	for (int ii = 0; ii < z_nodes.size(); ii++)
	{
		message += "Layer ID: " + std::to_string(ii) + "\n\n[ ";
		for (int jj = 0; jj < z_nodes[ii].size(); jj++)
			message += std::to_string(z_nodes[ii][jj]) + " ";
		message += "]\n-----------------------------------------\n";
	}
		
	message += "\n=============================================\n";

	if (print_to_terminal)
		std::cout << message << std::flush;
	if (out_file != NULL)
		(*out_file) << message << std::flush;

	return;

}


/*! \brief Function to print rho-node information for debugging and logging purposes.*/
void LayerManager::PrintNodeData_rho(std::ofstream *out_file, bool print_to_terminal)
{

	std::string message;

	message = "\n=============================================\n";
	message += "Node data (rho-nodes)\n";
	message += "=============================================\n\n";

	int N_nodes = rho_nodes.size();
	
	message += "Number of nodes: " + std::to_string(N_nodes) + "\n";

	message += "[ ";
	for (int ii = 0; ii < rho_nodes.size(); ii++)
		message += std::to_string(rho_nodes[ii]) + " ";
	
	message += "]\n=============================================\n";

	if (print_to_terminal)
		std::cout << message << std::flush;
	if (out_file != NULL)
		(*out_file) << message << std::flush;

	return;

}


// ==================================================================================
// Helpers - layer management
// ==================================================================================

/*! \brief Function to combine consecutive layers that have the same electrical properties within some tolerance.*/
void LayerManager::MergeLayersWithSameMaterial(double tol)
{

	// Initialize a vector to keep track of same-material layers
	
	std::vector<std::vector<int>> to_be_combined (1);
	to_be_combined.back().push_back(0);
	
	// Figure out which layers are to be merged
	
	for (int ii = 1; ii < layers.size(); ii++)
	{

		int jj = ii - 1; // The layer above

		// Check if the layer above is of the same material as this one
		if (std::abs(layers[ii].epsr - layers[jj].epsr) < tol &&
		    std::abs(layers[ii].mur - layers[jj].mur) < tol &&
		    std::abs(layers[ii].sigma - layers[jj].sigma) < tol &&
		    std::abs(layers[ii].sigmamu - layers[jj].sigmamu) < tol)
		{
			// It is the same, so remember that we need to combine this layer with the one above it
			to_be_combined.back().push_back(ii);
		}
		else
		{
			// It's not the same, so make a new subset of layers
			to_be_combined.push_back(std::vector<int> ());
			to_be_combined.back().push_back(ii);
		}

	}

	
	// DEBUGGING: Make sure the above works correctly
	// for (int ii = 0; ii < to_be_combined.size(); ii++)
	// {
	// 	std::cout << "Layer subset " << ii << " contains layers ( " << std::flush;
	// 	for (int jj = 0; jj < to_be_combined[ii].size(); jj++)
	// 		std::cout << to_be_combined[ii][jj] << " " << std::flush;
	// 	std::cout << ")\n" << std::endl;
	// }

	
	// Merge the layers
	std::vector<Layer> new_layers (to_be_combined.size());

	for (int ii = 0; ii < to_be_combined.size(); ii++)
	{

		// Copy over the first layer of this subset of the old layer-set
		new_layers[ii] = layers[to_be_combined[ii][0]];

		// The zmax of this combined layer is the same as the first layer of the subset. So are its material properties. But its zmin and height are different.
		new_layers[ii].zmin = layers[to_be_combined[ii].back()].zmin;
		new_layers[ii].h = new_layers[ii].zmax - new_layers[ii].zmin;

	}

	// Update the final list of layers
	layers.clear();
	layers = new_layers;
	new_layers.clear();


	// DEBUGGING: Make sure the above works correctly
	// for (int ii = 0; ii < layers.size(); ii++)
	// {
	// 	std::cout << "Layer " << layers[ii].LayerID << " from " << layers[ii].zmin << " to " << layers[ii].zmax << " with eps_r " << layers[ii].epsr << " and " << layers[ii].ObjectIDs.size() << " Objects." << std::endl;
	// }
	

	return;

}


/*! \brief Return the lower z-coordinate of the given layer. For the bottom-most layer, this function returns a coordinate a large distance below the stack-up. The value wouldn't matter in the context of the MGF computations.*/
double LayerManager::GetZmin(int idx)
{
	if (idx < -1 || idx > (int)layers.size() + 1)
		throw std::out_of_range("[ERROR] LayerManager::GetZmin(): Index out of bounds.");
	else if (idx == -1)
		return layers[0].zmax;
	else if (idx == (int)layers.size())
		return (layers.back().zmin - (layers[0].zmax - layers.back().zmin)*10.0);
	else
		return layers[idx].zmin;
}


/*! \brief Return the upper z-coordinate of the given layer. For the top-most layer, this function returns a coordinate a large distance above the stack-up. The value wouldn't matter in the context of the MGF computations.*/
double LayerManager::GetZmax(int idx)
{
	if (idx < -1 || idx > (int)layers.size() + 1)
		throw std::out_of_range("[ERROR] LayerManager::GetZmax(): Index out of bounds.");
	else if (idx == -1)
		return layers[0].zmax + (layers[0].zmax - layers.back().zmin)*10.0;
	else if (idx == (int)layers.size())
		return layers.back().zmin;
	else
		return layers[idx].zmax;
}


/*! \brief Return the height of the given layer.*/
double LayerManager::GetHeight(int idx)
{
	if (idx < -1 || idx > (int)layers.size() + 1)
		throw std::out_of_range("[ERROR] LayerManager::GetHeight(): Index out of bounds.");
	else if (idx == -1)
		return inf;
	else if (idx == (int)layers.size())
		return inf;
	else
		return layers[idx].h;
}






