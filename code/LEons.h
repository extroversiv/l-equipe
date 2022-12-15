/*
 * LEons.h
 *
 *  Created on: 09.01.2012
 *      Author: mik
 */

#ifndef LEONS_H_
#define LEONS_H_

#include <iostream>
#include <cmath>
#include <vector>
#include "reell.h"
#include "gsl/gsl_rng.h"
#ifdef BLAS
	#include "gsl/gsl_blas.h"
#endif
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;

class LE_ons {
	public:
		LE_ons(int, int, int, int);
		virtual ~LE_ons();

		int get_COL(){return COL;}
		int get_DIM(){return DIM;}
		int get_ROW(){return ROW;}
		reell get_normONS(int n){return normONS[n];}
		reell get_normCLV(int n){return normCLV[n];}

		void initONS(int seed);
		void orthonormalizeONS(int);

		bool saveProjections;

		unsigned get_numberSavedProjections(){return allSavedProjections.size();}

		void initCLV();
		void normalizeCLV(int);
		void backwardIterationCLV(unsigned, int);

		reell** ons;												//!< othonormal system for the Lypunov exponent calculation

	private:

		int myID, nP;

		gsl_rng* rng;												//!< random number generator

		int COL;													//!< number of columns (= vectors = number of Lyapunov exponents)
		int ROW;													//!< number of rows (= dimension = number of neurons*dimension of state variable)
		int DIM;													//!< the dimension of the neuron model's state variable (e.g. for phase neruons it's 1)

		vector<reell> normONS, projections;							//!< norm and projection vector used in orthogonalization

#ifdef PAR
		vector<reell> projectionsLocal;
#endif


		vector<vector<reell> > clv;									//!< coefficients of the covariant Lyapunov vectors w.r.t. to the orthonormal system
		vector<reell> normCLV;										//!< norms of covariant Lyapunov vecctors
		vector<vector<reell> > projStore;							//!< projection storage container for one orthonormalization
		vector<vector<vector<reell> > > allSavedProjections;		//!< array in which all projections are stored as long as the public saveProjections flag is set to true.
};

#endif /* LEONS_H_ */
