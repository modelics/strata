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

/************************ utility.hpp ************************

 * General-purpose utility functions.
 *
 * Author: Shashwat Sharma
 * Created on: June 19, 2019

 ***************************************************************/


#ifndef UTILITY_H
#define UTILITY_H

#include <complex>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unistd.h>
#include <chrono>

#if defined(_OPENMP)
	#include <omp.h>
#endif

namespace strata
{

	/*! \brief Function to emulate Matlab's linspace.*/
	template<typename T>
	inline void linspace(T start, T stop, int N, std::vector<T> &output)
	{

		output.resize(N);

		if (N == 1)
		{
			output[0] = start;
			return;
		}

		T dx = (stop - start)/((T)(N - 1));
		
		for (int ii = 0; ii < N; ii++)
			output[ii] = start + ((T)ii)*dx;

		return;

	}


	/*! \brief Function to emulate Matlab's logspace.*/
	template<typename T>
	inline void logspace(T start, T stop, int N, std::vector<T> &output)
	{

		linspace(start, stop, N, output);
		
		for (int ii = 0; ii < N; ii++)
			output[ii] = std::pow(10.0, output[ii]);

		return;

	}


	/*! \brief Function to normalize a vector so that its values lie between 0 and 1. Does not apply to complex numbers.*/
	template<typename T>
	inline void normalize_01(std::vector<T> &x)
	{

		// Find minimum element
		auto it_min = std::min_element(x.begin(), x.end());

		// Subtract smallest element
		std::for_each(x.begin(), x.end(), [it_min](T &d) { d -= *it_min;});

		// Find maximum element
		auto it_max = std::max_element(x.begin(), x.end());

		// Normalize to largest element
		std::for_each(x.begin(), x.end(), [it_max](T &d) { d /= *it_max;});

		return;

	}


	/*! \brief Function to normalize a vector to its maximum absolute value.*/
	template<typename T>
	inline void normalize_abs(std::vector<T> &x)
	{

		T max_val = 0.0;

		for (int ii = 0; ii < x.size(); ii++)
			if (std::abs(x[ii]) > std::abs(max_val))
				max_val = x[ii];

		for (int ii = 0; ii < x.size(); ii++)
			x[ii] /= max_val;

		return;

	}


	/*! \brief Function to reorder the elements of a given vector using a list of indices provided in another vector. Source: https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices*/
	template<class T>
	void reorder(std::vector<T> &v, std::vector<int> const &order)
	{   
		for (int s = 1, d; s < order.size(); ++s)
		{
			for (d = order[s]; d < s; d = order[d]);
			if (d == s)
				while (d = order[d], d != s)
					std::swap(v[s], v[d]);
		}

		return;
	}


	/*! \brief Function to print out the elements of an std::vector, assuming that the ofstream operator << is defined for the datatype stored in the vector.*/
	template<typename T>
	inline void PrintVector(std::vector<T> &x, std::string message = " ")
	{

		std::cout << "============ PrintVector() called ============" << std::endl;

		if (message != " ")
			std::cout << message << std::endl;
	
		for (int ii = 0; ii < x.size(); ii++)
			std::cout << x[ii] << std::endl;

		std::cout << "==============================================" << std::endl;
	
		return;

	}


	/*! \brief Function to print out the elements of an std::array, assuming that the ofstream operator << is defined for the datatype stored in the vector.*/
	template<typename T, size_t N>
	inline void PrintArray(std::array<T, N> &x, std::string message = " ")
	{

		std::cout << "============ PrintArray() called ============" << std::endl;

		if (message != " ")
			std::cout << message << std::endl;
	
		for (int ii = 0; ii < N; ii++)
			std::cout << x[ii] << std::endl;

		std::cout << "==============================================" << std::endl;
	
		return;

	}


	/*! \brief Function to print out the elements of a raw pointer array, assuming that the ofstream operator << is defined for the datatype stored in the array.*/
	template<typename T>
	inline void PrintArray(T *x, int N, std::string message = " ")
	{

		std::cout << "============ PrintArray() called ============" << std::endl;

		if (message != " ")
			std::cout << message << std::endl;

		for (int ii = 0; ii < N; ii++)
			std::cout << x[ii] << std::endl;

		std::cout << "==============================================" << std::endl;
	
		return;

	}


	/*! \brief Function to print out the keys and values of an std::map, assuming that the ofstream operator << is defined for the datatypes stored in the map.*/
	template<typename T1, typename T2>
	inline void PrintMap(std::map<T1, T2> &x, std::string message = " ")
	{

		std::cout << "============ PrintMap() called ============" << std::endl;

		if (message != " ")
			std::cout << message << std::endl;

		for (auto it = x.begin(); it != x.end(); it++)
			std::cout << "[" << it->first << ", " << it->second << "]" << std::endl;	

		std::cout << "===========================================" << std::endl;
	
		return;

	}


	/*! \brief Function that resembles std::to_string but allows the user to set the desired precision. Adapted from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values.*/
	template <typename T>
	inline std::string to_string_with_precision(const T a_value, const int n = 8)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << std::scientific << a_value;
		return out.str();
	}


	/*! \brief Function that resembles std::to_string but always uses scientific notation. Trailing zeros are removed.*/
	template <typename T>
	inline std::string to_string_scientific(const T a_value)
	{
		std::ostringstream out;
		out << std::scientific << a_value;
		return out.str();
	}

	
	/*! \brief Convenience function to convert a floating point number to a string without any decimals, for example to get (part of) a valid filename, with n digits of accuracy. E.g. if n = 3 and the number is 123456.789, the returned string will be 123E3.*/
	template <typename T>
	inline std::string to_string_without_decimal(const T a, const int n = 3)
	{
		std::stringstream stream;
		stream << std::fixed << std::setprecision(0) << a;
		std::string result = stream.str();
		if (result.length() > n)
			result = result.substr(0, n) + "E" + std::to_string((int)result.length() - n);
		return result;
	}


	/*! \brief The following two functions are adopted from https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process, in order to compute RAM usage by the program over the course of running it, for in-house profiling.*/
	inline int ParseLine(char* line)
	{
		// This assumes that a digit will be found and the line ends in " Kb".
		int i = std::strlen(line);
		const char* p = line;
		while (*p <'0' || *p > '9')
			p++;
		line[i-3] = '\0';
		i = atoi(p);
		return i;
	}


	/*! \brief Interface function to get RAM usage in kB.*/
	inline double GetMemoryUsage()
	{
		// Note: return value is in KB
	
		FILE* file = fopen("/proc/self/status", "r");
		int result = -1;
		char line[128];

		while (fgets(line, 128, file) != NULL)
		{
			if (strncmp(line, "VmRSS:", 6) == 0)
			{
				result = ParseLine(line);
				break;
			}
		}
		fclose(file);
	
		return (double) result;
	}


	/*! \brief Interface function to get the total system RAM in B.*/
	inline unsigned long long GetMemoryTotal()
	{
		// Note: return value is in B

		long pages = sysconf(_SC_PHYS_PAGES);
		long page_size = sysconf(_SC_PAGE_SIZE);
		return pages*page_size;
	}


	/*! \brief Interface function to get available RAM in kB.*/
	inline double GetMemoryAvailable()
	{
		// Note: return value is in kB
		return ((double)GetMemoryTotal()/1000.0) - GetMemoryUsage();
	}


	/*! \brief Class to control the display and logging of simulation information and progress, both for terminal and log files. Can also be used to log errors and warnings. Also logs profile data and can be used as a timer.*/
	class logstream
	{
	public:

		logstream(std::string filename_log = "./logfile.log", std::string filename_prof = "./profile.prof") { ResetLogFile(filename_log); ResetProfilingFile(filename_prof); return; };

		/// Reset log file name.
		void ResetLogFile(std::string const &filename) { print_to_log.open(filename); InitializeLog(); return; };

		/// Reset profiling file name.
		void ResetProfilingFile(std::string const &filename) { print_to_prof.open(filename); InitializeProfiling(); return; };
	
		/// Reset the list of warnings.
		void ResetWarnings() { warnings.clear(); return; };

		/// Reset the list of errors.
		void ResetErrors() { errors.clear(); return; };

		/// Log a new warning.
		void LogWarning(std::string message)
		{
			Log(message);
			warnings.push_back(message);
			return;
		};

		/// Log and display a new warning.
		void LogAndDisplayWarning(std::string message)
		{
			LogAndDisplay(message, true);
			warnings.push_back(message);
			return;
		};

		/// Log a new error.
		void LogError(std::string message)
		{
			Log(message);
			errors.push_back(message);
			return;
		};

		/// Log and display a new error.
		void LogAndDisplayError(std::string message)
		{
			LogAndDisplay(message, true);
			errors.push_back(message);
			return;
		};

	
		// ====== Functions to modify the log file ======

		/// Initialize the log.
		void InitializeLog()
		{
			if (!print_to_log)
			{
				std::cout << "[ERROR] Log file could not be opened or created." << std::endl;
				return;
			}

			auto now = std::chrono::system_clock::now();
			std::time_t time_now = std::chrono::system_clock::to_time_t(now);
			std::string time_str = std::ctime(&time_now);
			time_str.pop_back();

			LogAndDisplay("\n====================================================");
			LogAndDisplay("Began log on " + time_str);
			LogAndDisplay("====================================================\n");
		
			return;
		};

		/// Initialize the profiling file.
		void InitializeProfiling()
		{
			if (!print_to_prof)
			{
				std::cout << "[ERROR] Profiling file could not be opened or created." << std::endl;
				return;
			}

			auto now = std::chrono::system_clock::now();
			std::time_t time_now = std::chrono::system_clock::to_time_t(now);
			std::string time_str = std::ctime(&time_now);
			time_str.pop_back();

			print_to_prof << "\n====================================================" << std::endl;
			print_to_prof << "Began profiling on " << time_str << std::endl;
			print_to_prof << "====================================================\n" << std::endl;
			print_to_prof << "Entry" << "\t" << "Value" << "\t" << "Unit" << std::endl;
		
			return;
		};

		/// Print a message to terminal and to the log file.
		void LogAndDisplay(std::string message, bool print_to_terminal = true)
		{
			print_to_log << message << std::endl;
			if (print_to_terminal)
				std::cout << message << std::endl;
			return;
		};

		/// Print a message to terminal.
		void Display(std::string message)
		{
			std::cout << message << std::endl;
			return;
		};

		/// Print a message to terminal and to the log file.
		void Log(std::string message)
		{
			print_to_log << message << std::endl;
			return;
		};

		/// Print a message to terminal and to the log file.
		void AddProfileEntry(std::string entry, double value, std::string units)
		{
			print_to_prof << entry << "\t" << value << "\t" << units << std::endl;
			return;
		};
	

		/// \todo Destructor displays all errors and warnings that occurred, and also appends this information to the end of the log file.
		~logstream() {};
	

		std::vector<std::string> errors, warnings;
		std::ofstream print_to_log, print_to_prof;

		
		/// ====== Timer ======

		std::vector<std::clock_t> t_start;
		std::clock_t t_stop;
		std::vector<double> wt_start;
		double wt_stop;

		struct times { double dt = 0.0, dwt = 0.0; };

		/// Begin timer
		void StartTimer()
		{
			t_start.push_back(clock());
#if defined(_OPENMP)
			wt_start.push_back(omp_get_wtime());
#endif
		};

		/// Stop timer and get elapsed time in seconds
		times StopTimer()
		{
			if (t_start.size() == 0)
				throw std::logic_error("[ERROR] logstream::StopTimer(): Timer was never started.");

			times t;
			
			t_stop = clock();
			t.dt = (t_stop - t_start.back())/(double)CLOCKS_PER_SEC;
			t_start.pop_back();

#if defined(_OPENMP)
			wt_stop = omp_get_wtime();
			t.dwt = wt_stop - wt_start.back();
			wt_start.pop_back();
#endif

			return t;
		};

		/// Stop timer and get elapsed time in seconds
		void StopTimerAndPrintToLog(std::string message, bool print_to_terminal = true)
		{
			times t = StopTimer();
			AddProfileEntry(message, t.dt, "s");
			LogAndDisplay("[CPU TIME] " + message + ": " + to_string_with_precision(t.dt) + " s.", print_to_terminal);
#if defined(_OPENMP)
			LogAndDisplay("[WALL TIME] " + message + ": " + to_string_with_precision(t.dwt) + " s.", print_to_terminal);
#endif
		}

	};

}

#endif

