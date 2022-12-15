/*
 * structures.h
 *
 *  Created on: 31.01.2012
 *      Author: mik
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include "reell.h"

using namespace std;

struct st_spike
{
	reell time; 					// time until next spike (time interval)
	int neuron;						// spiking neuron
};

struct st_synapse{
	int post;						// postsynaptic neuron
	double cplg;					// coupling strength
	double prob;					// release probability
};

struct st_double_int
{
	double d;
	int i;
};

struct st_twoDlinear
{
	double alpha;
	double beta;
	double gamma;
	double delta;
	double Cw;
	double tauS;
};

struct st_in
{
	//The netcdf file stores doubles, therefore doubles need to be read in, even when reell is set to long double

	//! Either the number of spikes \a SC or the simulation time \a TC must be given. If both are given, then time \a TC is chosen.
	long long SR, SW, SC;											//!< spikes in rate finding(R), warmup(W) and calculation(C)
	double TR, TW, TC;												//!< time in rate finding(R), warmup(W) and calculation(C)

	double rateWnt, pR;												//!< wanted average firing rate in the network with precision pR

	int Nall;														//!< number of neurons
	int Nloc;														//!< number of local neurons on one node
	int Noffset;													//!< offset between local and global neuron index

	int homogNet;													//!< homogeneous network or not

	vector<int> neuronType;											//!< neuron types in network
	vector<double> rapidness;										//!< AP onset rapidness of the neurons
	vector<double> type1type2Para;									//!< parameter to interpolate btw. type 1 and type 2 neurons
	vector<double> reset;
	vector<double> threshold;
	vector<st_twoDlinear> twoDlinearParas;							//!< parameters for the twoDlinear model

	vector<double> tauM;											//!< membrane time constants of neurons
	vector<double> Iext;											//!< initial external current to the neurons
	vector<vector<reell> > init;									//!< initial states of the neurons
	vector<vector<st_synapse> > synapses;							//!< postsynaptic connections of the neurons

	vector<int> train;												//!< calculate the spike train of these neurons

	vector<int> ISIneurons;											//!< calculate the inter spike interval (ISI) statistics of these neurons
	int ISIstats;													//!< calculate this number of 'kinda' moments (1=mean, 2=cv, 3=skewness, 4=kurtosis)
	int ISIbins;													//!< number of bins for the ISI statistics distribution

	int LyapunovExponents;											//!< calculate this number of Lyapunov exponents
	int long long SWONS;											//!< spikes of the warmup of the orthonormal system before the Lyapunov exponents calculation
	int seedONS;													//!< seed for random number generator when creating the initial orthonormal system
	int ONstep;														//!< step size of orthonormalizations
	int LyapunovExponentsConvergence;								//!< Save the Lyapunov exponents at each orthonormalization step to a netcdf file.
	double pLEONS;														//!< precision of Lyapunov exponents

	int CLV;														//!< true false to calculate the covariant Lyapunov vectors
	int long long SWCLV;											//!< spikes of the warmup of the covariant Lyapunov vectors

	int saveFinalState;												//!< save the final state of the network

	int addCur;														//!< Add an additional time varying current? 0=false, 1=true
	int addCurHomo;													//!< Is this current the same for all neurons? 0=flase, 1=true
	vector<int> addCurNeurons;										//!< These neurons will receive the additional current.
	vector<double> addCurTimes;										//!< vector with the times at which the additional current is applied
	vector<vector<double> > addCurIext;								//!< vector with the additional currents

	int phases;														//! export the phases of all neurons?
	int distances;													//! export the distance between reference and perturbed trajectory (int value represents the norm)?

	int pertSpike;													//!< skip one spike
	int pertSynapse;												//!< skip one synaptic transmission
	double pertSize;												//!< perturbation size of given perturbation
	vector<double> pertVector;										//!< perturbation vector

	vector<double> phaseTimes;										//!< times at which the phases are saved
	vector<int> phaseSpikes;										//!< spikes at which the phases are saved

	double instPopRateBinSize;										//!< bin size with which to compute the instantaneous population firing rateC
	
};

struct st_out
{
	//Here reell should stay reell, because the variables might be used in the analysis.
	int N;															//!< number of neurons
	int long long SW, SC;											//!< spikes in warmup and calculation
	reell TW, TC;													//!< time in warmup and calculation													//! biological time of the simulation

	reell rateW, rateC;												//!< network-averaged firing rate in warmup and simulation

	vector<reell> finalCurrents;									//!< final external currents of the neurons
	vector<vector<reell> > finalStates;								//!< final states of the neurons

	vector<st_spike> spikeTrain;									//!< spike time and neuron, defined in LEnetwork.h

	int long long SWONS;											//!< spikes of the warmup of the orthonormal system before the Lyapunov exponents calculation
	int ONstep;														//!< step size of orthonormalizations
	vector<reell> LyapunovExponentsONS;								//!< Array of the asymptotic Lyapunov exponents.
	double pLEONS;													//!< precision of Lyapunov exponents
	vector<vector<reell> > LEconvergence;							//!< Save the Lyapunov exponents at each orthonormalization step,
	vector<reell> LEtimes;											//!< and save the times at which the exponents are saved.

	int long long SWCLV;											//!< spikes of the warmup of the covariant Lyapunov vectors
	vector<reell> LyapunovExponentsCLV;								//!< Array of the asymptotic Lyapunov exponents.
	vector<vector<reell> > localLyapunovExponents;							//!< Array of the local Lyapunov exponents at each orthonormalization (times stored in\a LEtimes).
	
	vector<reell> rateNeurons;										//!< the kinda moments of the inter spike interval distribution
	vector<reell> cvNeurons;
	vector<reell> skewnessNeurons;
	vector<reell> kurtosisNeurons;

	vector<vector<reell> > rateDist;								//!< and their distributions
	vector<vector<reell> > cvDist;
	vector<vector<reell> > skewnessDist;
	vector<vector<reell> > kurtosisDist;

	vector<vector<reell> > addCurRate;								//!< rates of individually chosen neurons

	int phases;														//! export the phases of all neurons?
	vector<reell> distances;										//! the distances between reference and perturbed trajectory (int value represents the norm)?

	vector<reell> phaseTimes;										//!< times at which the phases are saved
	vector<vector<reell> > phaseNeurons;							//!< phases of the neurons at phaseTimes
	
	vector<reell> instPopRate;							//!<time series of instaneous population firing rate
};

#endif /* STRUCTURES_H_ */
