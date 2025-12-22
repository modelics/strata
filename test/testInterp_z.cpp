#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <omp.h>
#include "MGF.hpp"


int main(int argc, char** argv)
{
    omp_set_num_threads(1);
    bool print_log = true;
    std::cout << "===========================" << std::endl;
    std::cout << "TestInterp()" << std::endl;
    std::cout << "===========================" << std::endl;

    std::string tech_file, out_file1, out_file2;

    if (argc < 2)
    {
        throw std::runtime_error("[ERROR] Layer file was not provided.");
    }
    else if (argc < 3)
    {
        tech_file = argv[1];
        out_file1 = "MGFdata_interp_real.txt";
        out_file2 = "MGFdata_interp_imag.txt";
    }
    else
    {
        tech_file = argv[1];
        out_file1 = argv[2];
        out_file2 = argv[3];
    }


    // ====== Stackup ======

    // Create a layer management object and parse the layer file
    LayerManager lm;
    lm.ProcessTechFile(tech_file);

    // Set the analysis frequency
    double f = 10e9;

    // Precompute frequency-dependent layer data
    lm.ProcessLayers(f);

    // Print layer data to the terminal for verification
    lm.PrintLayerData(NULL, true);


    // ====== Set up the source and observation points ======

    // Set the source and observation (obs) points

    double x_src = 0.0, y_src = 0.0;
    double x_obs = 5.0e-4, y_obs = 0.0;

    int Nz_interpolated = 20; // Number of points to be interpolated

    double dis_threshold = 1e-8;
    double z_min = lm.layers.back().zmin + dis_threshold;
    double z_max = lm.layers.front().zmax - dis_threshold;

    double z_interp_min = z_min;
    double z_interp_max = z_max;

    // We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
    std::vector<double> x_vec;
    std::vector<double> z_src_vec;
    std::vector<double> z_obs_vec;
    strata::linspace(z_interp_min, z_interp_max, Nz_interpolated, z_src_vec);
    strata::linspace(z_interp_min, z_interp_max, Nz_interpolated, z_obs_vec);


    // ====== Create the interpolation grid ======

    // The interpolation grid consists of nodes along the z and rho axes.
    // When computing the MGF for arbitrary source and observation coordinates:
    // - along the z axis, Lagrange polynomial interpolation of a given order is performed,
    // - along the rho axis, Lagrange polynomial interpolation of a given order is performed.

    // In this example, Nz nodes are used to create the interpolation table
    // ====== Initialize the MGF class ======

    // This class stores all the settings we want to use
    MGF_settings s;
    s.interpolate_z = true;
    s.update_z_nodes = true;

    // Strata comes with Matlab-like linspace and logspace functions for convenience
    std::vector<double> rho_nodes = {0.0, 5.11155e-05, 0.000102231, 0.000153347, 0.000204462, 0.000255578, 0.000306693, 0.000357809, 0.000408924, 0.00046004, 0.000511155, 0.000562271, 0.000613386, 0.000664502, 0.000715617, 0.000766733, 0.000817849, 0.000868964, 0.00092008, 0.000971195};
    //std::vector<double> rho_nodes = {0.0, 2.5e-3, 5.0e-3};

    // Provide the layer manager with the z and rho nodes
    lm.z_nodes.resize(lm.layers.size());
    lm.z_nodes[0] = {0.000078,0.000083,0.000088};
    lm.z_nodes[1] = {0.000071,0.000074,0.000078};
    lm.z_nodes[2] = {0.000064,0.000067,0.000071};
    lm.z_nodes[3] = {0.000057,0.000060,0.000064};
    lm.z_nodes[4] = {0.000040,0.000048,0.000057};
//    lm.z_nodes[0] = {0.000078,0.000088};
//    lm.z_nodes[1] = {0.000071,0.000078};
//    lm.z_nodes[2] = {0.000064,0.000071};
//    lm.z_nodes[3] = {0.000057,0.000064};
//    lm.z_nodes[4] = {0.000040,0.000057};


    lm.ClearNodes_rho();
    lm.InsertNodes_rho(rho_nodes);

    // This class is the "engine" which will compute the MGF
    MGF mgf;

    // Tell the class that we want to use the interpolation method.
    s.method = MGF_INTERPOLATE;
    s.filename = "./interpolation_table.txt";
    s.export_table = false;

    // For source and observation points close to each other, the interpolation can be inaccurate due to the singularity of the MGF.
    // So, instruct the MGF engine to extract the quasistatic part in the spectral domain, and then add it back in the spatial domain analytically.
    s.extract_quasistatic = true;

    // Polynomial interpolation order
    s.order = 2;
    s.order_z = 2;

    // Initialize the class with the given stackup and chosen settings.
    // The interpolation table will be generated at this point, so the initialization step could take a while.
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    mgf.Initialize(f, lm, s);
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    static int Interpolation_initialization = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count();

    // ====== Compute the MGF ======

    // Create an output file where the MGF will be exported for post-processing
    std::ofstream outfile1(out_file1);
    std::ofstream outfile2(out_file2);

    // For post-processing, we'll store the frequency and positions along the z axis in the header
    outfile1 << "Frequency: " << f << " Hz" << std::endl;
    outfile2 << "Frequency: " << f << " Hz" << std::endl;

    // The lateral separation between the source and observation points (rho) and all the MGF components is tabulated for each observation point
    outfile1 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;
    outfile2 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

    std::cout << "Number of rho points: " << rho_nodes.size() << std::endl;

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();
    mgf.AdaptiveInterpolation();
    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    static int AdaptiveInterpolation = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count();

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    for (int ii = 0; ii < Nz_interpolated; ii++) {
        for (int jj = 0; jj < Nz_interpolated; jj++) {

            // The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
            // The working source and observation layers can be found with the FindLayer() method of the layer manager.
            // In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
            int i = lm.FindLayer(z_src_vec[ii]);
            int m = lm.FindLayer(z_obs_vec[jj]);
            mgf.SetLayers(i, m); // Source first, observation second


            // In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
            double x_diff = x_obs - x_src;
            double y_diff = y_obs - y_src;

            // The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
            // In addition, there is a scalar component which is just a complex number.
            std::array<std::complex<double>, 9> G_dyadic;
            std::complex<double> G_phi;

            // Compute the MGF
            mgf.ComputeMGF(x_diff, y_diff, z_obs_vec[jj], z_src_vec[ii], G_dyadic, G_phi);

            // The dyadic MGF components are now stored in G_dyadic, while the scalar MGF is stored in G_phi.


            // ====== Export data to the output text file ======

            // Lateral separation
            outfile1 << z_src_vec[ii] << " " << z_obs_vec[jj] << " " <<
                     std::real(G_dyadic[0]) << " " << std::real(G_dyadic[1]) << " " << std::real(G_dyadic[2]) << " " <<
                     std::real(G_dyadic[3]) << " " << std::real(G_dyadic[4]) << " " << std::real(G_dyadic[5]) << " " <<
                     std::real(G_dyadic[6]) << " " << std::real(G_dyadic[7]) << " " << std::real(G_dyadic[8]) << " "
                     << std::real(G_phi) << std::endl;

            outfile2 << z_src_vec[ii] << " " << z_obs_vec[jj] << " " <<
                     std::imag(G_dyadic[0]) << " " << std::imag(G_dyadic[1]) << " " << std::imag(G_dyadic[2]) << " " <<
                     std::imag(G_dyadic[3]) << " " << std::imag(G_dyadic[4]) << " " << std::imag(G_dyadic[5]) << " " <<
                     std::imag(G_dyadic[6]) << " " << std::imag(G_dyadic[7]) << " " << std::imag(G_dyadic[8]) << " "
                     << std::imag(G_phi) << std::endl;

        }
    }
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    static int MGF_computation = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count();

    if (print_log)
    {
        std::cout << "Nz: " << mgf.lm.z_nodes[0].size() << std::endl;
        std::cout << "Nz_interpolated: " << Nz_interpolated << std::endl;
        std::cout << "z_interp_min: " << z_interp_min << ", z_interp_max: " << z_interp_max << std::endl;
        std::cout << "Initialize the interpolation table: " << Interpolation_initialization << "(ms)" << std::endl;
        std::cout << "Update the interpolation table: " << AdaptiveInterpolation << "(ms)" << std::endl;
        std::cout << "Compute the MGF: " << MGF_computation << "(ms)" << std::endl;
    }


    return 0;

}






