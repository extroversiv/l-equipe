/*
 * LEneuron.h
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#ifndef LENEURON_H_
#define LENEURON_H_

#include "reell.h"
#include <iostream>
#include <vector>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_math.h"

using namespace std;

//! The neuron class.
/*!
  The base class of all neuron models.
*/

class LE_neuron
{
	public:
		LE_neuron(reell, vector<reell>, vector<reell>);
		virtual ~LE_neuron();

		reell get_spikeTime();
		int get_stateDim();
		vector<reell> get_state();
		virtual reell get_phase() = 0;
		reell get_externalCurrent();

		void set_state(vector<reell>&);
		virtual void set_phase(reell&) = 0;
		virtual void set_externalCurrent(reell&) = 0;


		virtual void dummy_evolve(reell, vector<reell>*);
		virtual void evolve_dt(reell, vector<reell>*) = 0;
		virtual void evolve_spike(reell) = 0;
		virtual void reset() = 0;
		virtual void dummy_PostUpdate();

		virtual vector<vector<reell> > calc_JacElem_postsynaptic(vector<reell>*, reell) = 0;
		virtual vector<vector<reell> > calc_JacElem_self(vector<reell>*) = 0;
		virtual vector<vector<reell> > calc_JacElem_spiking(vector<reell>*) = 0;
		virtual void dummy_Jacobian_spiking(vector<reell>*) = 0;			//!< Model dependent function to get a value (e.g. the velocity) of the spiking neuron needed for the single spike Jacobian.

	protected:
		virtual void calc_spikeTime() = 0;

		vector<reell> state;											//!< the state variable can be multidimensional in derived classes -> pointer
		vector<reell> stateReset;										//!< reset value
		vector<reell> stateThreshold;									//!< threshold value

		reell tauM;														//!< the membrane time constant
		reell spikeTime;												//!< the next spike time of the neuron
		reell Irheo;													//!< the rheo base current (the minimum injected current for the neuron to start tonic firing)
		reell Iext;														//!< the external current

		bool calcSpikeTime;
};

#endif /* LENEURON_H_ */
