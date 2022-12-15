/*
 * LEneuron.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#include "LEneuron.h"

LE_neuron::LE_neuron(reell tau, vector<reell> reset, vector<reell> threshold)
{
	//! All time related quantities are expressed in terms of tauM internally.
	tauM = tau;

	//! Every neuron has some reset and threshold values.
	stateReset = reset;
	stateThreshold = threshold;

	Iext = 0;

	calcSpikeTime = false;
}

LE_neuron::~LE_neuron()
{

}

reell LE_neuron::get_externalCurrent()
{
	return Iext;
}


reell LE_neuron::get_spikeTime()
{
	//! Returns the next spike time of the neuron.
	//! If it is not up to date, then it is calculated first.
	if (calcSpikeTime)
		calc_spikeTime();

	if (spikeTime < 0)
	{
		cout << "The neurons next spike time " << spikeTime << " was in the past!" << endl;
		throw(1);
	}

	return spikeTime*tauM;
}

int LE_neuron::get_stateDim()
{
	return state.size();
}

vector<reell> LE_neuron::get_state()
{
	return state;
}

void LE_neuron::set_state(vector<reell>& newState)
{
	state = newState;
	calcSpikeTime = true;
}


//! This function can be overloaded by the derived classes to store some dummy variables.
void LE_neuron::dummy_evolve(reell dt, vector<reell>* dummy)
{

}

void LE_neuron::dummy_PostUpdate()
{

}

