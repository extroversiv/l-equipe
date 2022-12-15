/*
 * LEnetwork.h
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#ifndef LENETWORK_H_
#define LENETWORK_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include "gsl/gsl_rng.h"
#include "reell.h"
#include "structures.h"
#include "LEneuron.h"
#include "LErapidtheta.h"
#include "LElif.h"
#include "LEtype1type2.h"
#include "LEclif.h"
#include "LEtwoDlinear.h"
#include "LEons.h"
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;

//! The network class.
/*!
  The class that handles the neuronal network setup and evolution.
*/

class LE_network{
	public:
		LE_network(st_in*);
		virtual ~LE_network();

		void set_rngSynapses(unsigned);
		void set_externalCurrents_Nloc(st_in*);
		void set_externalCurrents_Nloc(vector<reell>&);
		void set_synapses(st_in*);
		void set_state_Nloc(vector<vector<reell> >&);
		void set_phase_Nloc(vector<reell>&);
		void set_externalCurrents_neurons(vector<int>&, vector<reell>&);

		vector<st_spike> get_spikes();
		int get_N();
		int get_Nloc();
		bool get_allNeuronSame();
		int get_stateDim(int);
		vector<vector<reell> > get_state_Nloc();
		vector<reell> get_phase_Nloc();
		vector<reell> get_externalCurrents_Nloc();
		vector<reell> get_externalCurrents_neurons(vector<int>&);

		reell find_nextSpikes();
		void evolve_dt(reell);
		void evolve_spike(bool);

		void multiply_JacobianONS(LE_ons*, int);

	private:
		inline void calcJacobian_and_updateNeurons();					//!< Calculates the Jacobian elements and updates the neurons.
		inline void updateNeurons();									//!< Updates the neurons when receiving a spike.

		//inline void init_randomgraph(reell, reell);					//!< initialize random graph (for test)

		inline int global2local(int);									//! transform from global neuron ID to the local neuron ID
		inline int local2global(int);									//! transform from local neuron ID to the global neuron ID
		inline int spikingNode(int);

		int myID, nP;

		int Nall, Nloc, Noffset;										//!< global number of all neurons, the number of the local neurons and the offset between local and global ID

		vector<vector<vector<vector<reell> > > > jacobian;				//!< storage container for the nonzero elements of the single spike Jacobians (multiple if sync spikes)
																		//!< outermost has dynamic length of the number of postsynaptic neurons + 1 for the spiking neuron
		vector<LE_neuron*> neurons;										//!< array of all neurons (possibly different models)

		vector<vector<st_synapse> >* synapses;							//!< connection matrix with postynaptic neurons and coupling strength

		vector<st_spike> spikes;										//!< vector with spikes in each interval

		bool allNeuronsSame;											//!< True if all neurons are of the same type (important for Jacobian)
		bool allPhaseNeurons;											//!< True if all neurons are phase neurons
		bool allClifNeurons;											//!< True if all neurons are cLIF neurons


		gsl_rng* rngSyn;												//!< random number generator for stochastic synapses
		vector<reell> pSyn;												//!< release probability of the synapses for each neuron


		// some dummy variables for the evolution and the single spike Jacobian
		vector<reell> dummyEvolve;
		vector<reell> dummyJacobian;
};

#endif /* LENETWORK_H_ */
