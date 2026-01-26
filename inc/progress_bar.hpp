/****************************** progress_bar.hpp *****************************

 * Class to allow printing progress bars to terminal.
 * With help from: https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
 *
 * Author: Shashwat Sharma
 * Created on: Nov 28, 2018

 *****************************************************************************/

/*

Usage:

// Initialize:
progress_bar bar (<N_iters>, <width of bar in # characters>, <message to display at the start>);

// Right after beginning the loop:
bar.UpdateProgressBar(<current loop index>);

// After the end of the loop:
bar.EndProgressBar(<message to display at the end>);

*/


#ifndef PROGBAR_H
#define PROGBAR_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>


/*! \brief Class to store and manipulate points or position vectors in 3D, in Cartesian coordinates.*/
class progress_bar
{
public:

	double progress = 0.0;
    double step = 0.0;
    int step_iter = 0;
    int N_iters = 0;
    int bar_width = 50;
	bool verbose = true;

	/*! \brief Constructor.*/
	progress_bar(int _N_iters, int _bar_width = 50, std::string message = "", bool _verbose = true)
    {
		verbose = _verbose;

		if (!verbose)
			return;

		N_iters = _N_iters;
        bar_width = _bar_width;
        progress = 0.0;

        // Steps in which progress is to be incremented
        step = 100.0/((double)bar_width);

        // Iteration steps at which to update progress
        step_iter = N_iters/bar_width;

		if (step_iter == 0)
		{
			step_iter = 1;
			step = 100.0/((double)N_iters) + 1.0;
		}

		std::cout << message << std::endl;

        return;
    };

    /*! \brief Function to update the progress bar at each iteration.*/
    void UpdateProgressBar(int iter)
    {

		if (!verbose)
			return;

        std::cout << "[";

        int pos = bar_width*(progress/100.0);

		if (iter + 1 == N_iters)
		{
			progress = 100.0;
			pos = bar_width;
		}

        for (int ii = 0; ii < bar_width; ++ii)
        {
            if (ii <= pos)
                std::cout << "=";
            // else if (ii == pos)
            //     std::cout << ">";
            else
                std::cout << " ";
        }
		
        std::cout << "] " << int(progress) << " %\r" << std::flush;

		
        if ((step_iter == 1) || (iter % step_iter == 0 && iter != 0 && progress < 100.0))
		{
            progress += step;
			progress = std::min(100.0, progress);
		}
		
		
        return;
		
    };

    void EndProgressBar(std::string message = "")
    {
		if (!verbose)
			return;

        std::cout << std::endl << message << std::endl;

		return;
    };

};


#endif
