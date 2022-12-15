/*
 * LEinputoutput.h
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */


#ifndef LEINPUTOUTPUT_H_
#define LEINPUTOUTPUT_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <netcdfcpp.h>
#include "structures.h"
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;

class LE_input_output
{
	public:
		LE_input_output();
		virtual ~LE_input_output();

		void read_netcdf(string, string, string, st_in*);
		void write_netcdf(string, st_out*, string, string, string);

		void display_in(st_in*);
		void display_out(st_out*);

		void gather_Nall(vector<double>*);
		void gather_Nall(vector<int>*);

	protected:
		int myID, nP;
};



#endif /* LEINPUTOUTPUT_H_ */
