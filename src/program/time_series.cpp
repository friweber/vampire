//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
//

// Standard Libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 


// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void time_series(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::time_series has been called" << std::endl;

	double temp=sim::temperature;

   // Set equilibration temperature only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
	   // Set equilibration temperature
	   sim::temperature=sim::Teq;
   }

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();
	}

   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set simulation temperature
	   sim::temperature = temp;

      // Reset mean magnetisation counters
      stats::mag_m_reset();

   }

   	// Vector to hold temperature data from file
	std::vector<double> temperature_data;
	std::ifstream file(std::getenv("TEMPERATURE_FILE")); // File containing temperature data
	double temp_list;

	// Read the temperature data from file
	while (file >> temp_list) {
		temperature_data.push_back(temp_list);
		
	}
	file.close();

	// Ensure there is temperature data
	if (temperature_data.empty()) {
		std::cerr << "Error: No temperature data loaded." << std::endl;
		exit(EXIT_FAILURE); // Exit the program if no data is loaded
	}

	size_t temp_index = 0; // Index to track the current temperature data point

	// Perform Time Series
	while (sim::time < sim::equilibration_time + sim::total_time) {

		// Update the system temperature from file data
		if (temp_index < temperature_data.size()) {
			sim::temperature = temperature_data[temp_index++];
		} else {
			std::cerr << "Warning: Not enough temperature data points; reusing the last available temperature." << std::endl;
			// Optional: Repeat the last temperature or implement other logic
			sim::temperature = temperature_data.back();
		}

		// Integrate system
		sim::integrate(1);

		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();
	}


}

}//end of namespace program
