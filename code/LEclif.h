/*
 * LEclif.h
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#ifndef LECLIF_H_
#define LECLIF_H_

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include "LEneuron.h"

using namespace std;

//! The correlated leaky integrate and fire class.

class LE_clif : public LE_neuron
{
	public:
		LE_clif(reell);
		virtual ~LE_clif();

		virtual reell get_phase();

		virtual void set_phase(reell&);
		virtual void set_externalCurrent(reell&);

		virtual void evolve_dt(reell, vector<reell>*);
		virtual void evolve_spike(reell);
		virtual void reset();
		virtual void dummy_PostUpdate();

		virtual vector<vector<reell> > calc_JacElem_postsynaptic(vector<reell>*, reell);
		virtual vector<vector<reell> > calc_JacElem_self(vector<reell>*);
		virtual vector<vector<reell> > calc_JacElem_spiking(vector<reell>*);
		virtual void dummy_Jacobian_spiking(vector<reell>*);

	protected:
		virtual void calc_spikeTime() = 0;
		virtual reell powGamma(reell) = 0;

		reell gamma;												//!< the ratio of membrane time constant and synaptic time constant (gamma = tauM/tauI)
		reell tauS;													//!< the synaptic time constant tauS = tauM/gamma;

		//! Sums up the time between two successive spikes in the network used for the jacobian calculation.
		//! Here, this is done for the exponential with the membrane time constant and with the synaptic time constant.
		reell jacISIexpM;
		reell jacISIexpS;

};


//!-------------------------- gamma = 2 --------------------------


class LE_clif_fast2 : public LE_clif
{
	public:
		LE_clif_fast2(reell);
		virtual ~LE_clif_fast2();

	protected:
		inline virtual void calc_spikeTime();
		inline virtual reell powGamma(reell);
};


//!------------------------ gamma = 1/2 --------------------------


class LE_clif_slow2 : public LE_clif
{
	public:
		LE_clif_slow2(reell);
		virtual ~LE_clif_slow2();

	protected:
		inline virtual void calc_spikeTime();
		inline virtual reell powGamma(reell);
};


//!-------------------------- gamma = 3 --------------------------


class LE_clif_fast3 : public LE_clif
{
	public:
		LE_clif_fast3(reell);
		virtual ~LE_clif_fast3();

	protected:
		inline virtual void calc_spikeTime();
		inline virtual reell powGamma(reell);

		reell MACHINE_PRECISION;
};


//!------------------------- gamma = 1/3 -------------------------


class LE_clif_slow3 : public LE_clif
{
	public:
		LE_clif_slow3(reell);
		virtual ~LE_clif_slow3();

	protected:
		inline virtual void calc_spikeTime();
		inline virtual reell powGamma(reell);

		reell MACHINE_PRECISION;
};


#endif /* LECLIF_H_ */
