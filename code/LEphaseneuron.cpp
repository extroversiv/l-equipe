/*
 * LEphaseneuron.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#include "LEphaseneuron.h"

LE_phaseneuron::LE_phaseneuron(reell tauM, vector<reell> reset, vector<reell> threshold) : LE_neuron(tauM, reset, threshold)
{
	//! The state variable is one-dimensional in the phase neuron class.
	state = vector<reell> (1);
}

LE_phaseneuron::~LE_phaseneuron()
{

}

reell LE_phaseneuron::get_phase()
{
	//! In the phase neuron class, the phase is obviously the state variable. This is returned here.
	return state[0];
}

void LE_phaseneuron::set_phase(reell& phase)
{
	//! If the phase is above threshold, it's set to the threshold.
	state[0] = (phase < stateThreshold[0]) ? phase : stateThreshold[0];

	calcSpikeTime = true;
}

inline void LE_phaseneuron::calc_spikeTime()
{
	spikeTime = (stateThreshold[0] - state[0])/w;

	calcSpikeTime = false;
}

//! The phase evolves with constant phase velocity w between spike events.
void LE_phaseneuron::evolve_dt(reell dt, vector<reell>* dummy)
{
	dt /= tauM;

	state[0] += w*dt;

	calcSpikeTime = true;
	//! We cannot just use spikeTime -= dt because this acccumulates precision losses. For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.
}

//! This is the neuron model specific phase transition curve.
void LE_phaseneuron::evolve_spike(reell c)
{
	state[0] = PTC(c);

	calcSpikeTime = true;
}


//! Set the state variable to the reset value.
void LE_phaseneuron::reset()
{
	state = stateReset;

	calcSpikeTime = true;
}

//! Calculate the elements using the derivative of the phase transition curve of the neuron model.
//! \awSpike is one-dimensional and contains the phase velocity of the spiking neuron
vector<vector<reell> > LE_phaseneuron::calc_JacElem_postsynaptic(vector<reell>* dummy, reell c)
{
	vector<vector<reell> > jacTmp(2, vector<reell>(1));				//state.size = 1 for phase neurons

	//! Two elements for the diagonal and the nondiagonal part, respectively. Each is one-dimensional for phase neurons.
	jacTmp[0][0] = PTCprime(c);										//!< the diagonal element
	jacTmp[1][0] = w/tauM/(*dummy)[0]*(1 - jacTmp[0][0]);					//!< the offdiagonal element

	return jacTmp;
}

vector<vector<reell> > LE_phaseneuron::calc_JacElem_self(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(1));

	jacTmp[0][0] = 1;												//!< the diagonal element

	return jacTmp;
}

vector<vector<reell> > LE_phaseneuron::calc_JacElem_spiking(vector<reell>* dummy)
{
	//! nothing special for the spiking neuron in phase neuron models,
	//! because the phase velocity is the same before and after the spike
	return calc_JacElem_self(dummy);
}

//! For phase neurons the phase velocity of the spiking neuron is needed to set up the Jacobian.
void LE_phaseneuron::dummy_Jacobian_spiking(vector<reell>* dummy)
{
	dummy->resize(1);

	(*dummy)[0] = w/tauM;
}
