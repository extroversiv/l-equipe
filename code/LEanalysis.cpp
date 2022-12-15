/*
 * LEanalysis.cpp
 *
 *  Created on: 04.01.2012
 *      Author: mik
 */

#include "LEanalysis.h"

LE_analysis::LE_analysis(LE_network* network) : net(network)
{
	//! initialize parallel environment
	nP = 1;																	//!< number of processors involved
	myID = 0;																//!< local ID of the processors, for unique communications

#ifdef PAR
	//! Initialize the parallel environment.
	MPI_Comm_size(MPI_COMM_WORLD,&nP);
	MPI_Comm_rank(MPI_COMM_WORLD,&myID);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	ons = new LE_ons(1, 1, 1, 1);
}

LE_analysis::~LE_analysis()
{
	delete ons;
}



void LE_analysis::multitask(st_in* in, st_out* out, bool applyPerturbation)
{
	//! This method runs the simulation and analysis. The different results requested are all calculated from the same trajectory.
	//! An time-varying external current can be applied, although in this case, the Lyapunov spectrum cannot be calculated.


	if (myID == 0)
	{
		cout << endl << "on the " << ((!applyPerturbation) ? "reference" : "perturbed") << " trajectory, gonna calculate: " << endl;

		if (in->phaseTimes.size()) cout << "\t * the phases at the provided times" << endl;
		if (in->phaseSpikes.size()) cout << "\t * the phases at the provided spikes" << endl;
		if (in->LyapunovExponents > 0) cout << "\t * " << in->LyapunovExponents << " Lyapunov exponents" << endl;
		if (in->CLV > 0) cout << "\t * " << in->LyapunovExponents << " covariant Lyapunov vectors" << endl;
		if (in->train.size() > 0) cout << "\t * the spike train of " << in->train.size() << " neurons" << endl;
		if ((in->ISIneurons.size() > 0) && (in->ISIstats > 0)) cout << "\t * the inter spike interval statistics of " << in->ISIneurons.size() << " neurons, the first " << in->ISIstats << " kinda moments" << endl;
	}


	this->applyPerturbation = applyPerturbation;


	//! run preprocessing first
	preprocessing(in, out);




	if (myID == 0) cout << endl << "simulation with SC = " << in->SC << " or TC = " << in->TC << " ms ..." << endl;

	out->TC = 0;
	out->SC = 0;

	int Tdone = 1, Sdone = 1;
	notConverged = false;
	while ((out->TC < in->TC) || (out->SC < in->SC) || notConverged)
	{
		prespike(in, out);

		updatespike(in, out);

		//! Display how much is done already. If longer simulations are needed too reach the desired precision, it will exceed 100%.
		if ((in->TC > 0) && (out->TC >= in->TC/10.*Tdone))
			if (myID == 0) cout << "\t" << Tdone++*10 << "% of TC done ... " << out->SC << " spikes @ t = " << out->TC/1000 << "s" << endl;

		if ((in->SC > 0) && (out->SC >= in->SC/10.*Sdone))
			if (myID == 0) cout << "\t" << Sdone++*10 << "% of SC done ... " << out->SC << " spikes @ t = " << out->TC/1000 << "s" << endl;
	}

	out->rateC = out->SC/out->TC/net->get_N();
	if (myID == 0) cout << "\tspikes: " << out->SC << "\ttime: " << out->TC/1000 << "s\t -> avg. rate: " << out->rateC*1000 << " Hz" << endl;



	//!run postprocessing last!
	postprocessing(in, out);


}



void LE_analysis::preprocessing(st_in* in, st_out* out)
{

	if (myID == 0) cout << endl << "preprocessing ... " << endl << endl;

	//***************
	//! check  stuff
	//***************

	//! The analysis runs until the maximum of (spikes\a SC or time\a TC is reached)
	if ((in->TC <= 0) && (in->SC <= 0))
	{
		if (myID == 0) cout << "SC = " << in->SC << "\tTC = " << in->TC << endl;
		if (myID == 0) cout << "Neither the number of spikes SC or the simulation time TC were provided! -->> exit" << endl;
		throw(1);
	}

	//! Check, whether the Jacobian needs to be calculated, e.g. for the calculation of the Lyapunov exponents
	calcJacobian = false;
	if (in->LyapunovExponents > 0) //&& others
		calcJacobian = true;

	if (calcJacobian)
	{
		if (!net->get_allNeuronSame())
		{
			if (myID == 0) cout << "Cannot calculate the single spike Jacobian because the neurons in the network are of different types." << endl;
			throw(1);
		}

		if (in->addCur)
		{
			if (myID == 0) cout << "If the external current is changed within a spike interval, what does that mean for the Jacobian?" << endl;
			throw(1);
		}
	}

	//********************
	//! initialize  stuff
	//********************

	//! Pointer to the orhonormal system, e.g. for the calculation of the Lyapunov exponents.
	if (in->LyapunovExponents > 0)
		pre_LyapunovExponents(in, out);

	if (in->ISIstats > 0)
		pre_ISIstatistics(in, out);

	//! This only works with the phase models so far.
	if (in->pertSize || in->pertSpike || in->pertSynapse)
		pre_perturbation(in, out);

	//! get the current external currents of the neurons that receive an additional current
	if (in->addCur)
	{
		curCount = 0;

		currentsOrg.resize(in->addCurNeurons.size());
		currentsTmp.resize(in->addCurNeurons.size());

		currentsOrg = net->get_externalCurrents_neurons(in->addCurNeurons);
	}

	//! initialize population firing rate times series vector
	if (in->instPopRateBinSize > 0)
		out->instPopRate = vector<double> (int(ceil(in->TC/in->instPopRateBinSize)));

}



void LE_analysis::prespike(st_in* in, st_out* out)
{
	//! if a perturbation was or will be applied (e.g. for distance measurements), save the phase of the neurons at the given times
	//! if an additional current is provided, change the external currents
	if (in->phaseTimes.size() || in->addCur)
	{
		//! The precision seems to be a problem, since there are a couple of difference calculations that influence the network evolution.
		reell timeToNextSpike = net->find_nextSpikes();
		reell timeToNextPhaseMeasurement = (in->phaseTimes.size() && (phaseTimeCount < in->phaseTimes.size())) ? in->phaseTimes[phaseTimeCount] - out->TC : timeToNextSpike + 1;
		reell timeToNextCurrentStep = (in->addCur && (curCount < in->addCurTimes.size())) ? in->addCurTimes[curCount] - out->TC : timeToNextSpike + 1;

		while((timeToNextSpike > timeToNextPhaseMeasurement) || (timeToNextSpike > timeToNextCurrentStep))
		{
			//! if the next spike time is after the next time to measure the phase, then evolve the network to the next \aphaseTime.
			//! if the next spike time is after the next current change, then evolve to the next \aaddCurTime.

			if (timeToNextPhaseMeasurement <= timeToNextCurrentStep)
			{
				if (timeToNextSpike > timeToNextPhaseMeasurement)
				{
					//! evolve the network until phaseTime
					reell dt = timeToNextPhaseMeasurement;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;

					out->phaseNeurons[phaseTimeCount] = net->get_phase_Nloc();
					if (myID == 0) out->phaseTimes[phaseTimeCount] = out->TC;

					//! go to the next time step
					phaseTimeCount++;
				}

				if (timeToNextSpike > timeToNextCurrentStep)
				{
					//! evolve the neurons to the time at which the current is changed
					reell dt = timeToNextCurrentStep - timeToNextPhaseMeasurement;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;

					//! change the current of the specified neurons
					for (unsigned n=0; n<in->addCurNeurons.size(); n++)
					{
						int nHomo = (in->addCurHomo) ? 0 : n;
						currentsTmp[n] = currentsOrg[n] + in->addCurIext[curCount][nHomo];
					}
					net->set_externalCurrents_neurons(in->addCurNeurons, currentsTmp);

					//! go to the next time step
					curCount++;
				}
			}
			else	//the order is the other way around
			{
				if (timeToNextSpike > timeToNextCurrentStep)
				{
					//! evolve the neurons to the time at which the current is changed
					reell dt = timeToNextCurrentStep;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;

					//! change the current of the specified neurons
					for (unsigned n=0; n<in->addCurNeurons.size(); n++)
					{
						int nHomo = (in->addCurHomo) ? 0 : n;
						currentsTmp[n] = currentsOrg[n] + in->addCurIext[curCount][nHomo];
					}
					net->set_externalCurrents_neurons(in->addCurNeurons, currentsTmp);

					//! go to the next time step
					curCount++;
				}

				if (timeToNextSpike > timeToNextPhaseMeasurement)
				{
					//! evolve the network until phaseTime
					reell dt = timeToNextPhaseMeasurement - timeToNextCurrentStep;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;

					out->phaseNeurons[phaseTimeCount] = net->get_phase_Nloc();
					if (myID == 0) out->phaseTimes[phaseTimeCount] = out->TC;

					//! go to the next time step
					phaseTimeCount++;
				}
			}

			timeToNextPhaseMeasurement = (in->phaseTimes.size() && (phaseTimeCount < in->phaseTimes.size())) ? in->phaseTimes[phaseTimeCount] - out->TC : timeToNextSpike + 1;
			timeToNextCurrentStep = (in->addCur && (curCount < in->addCurTimes.size())) ? in->addCurTimes[curCount] - out->TC : timeToNextSpike + 1;

		}
	}

	//! The normal evolution of all neurons without any external changes

	reell timeToNextSpike = net->find_nextSpikes();

	net->evolve_dt(timeToNextSpike);

	out->TC += timeToNextSpike;

}

void LE_analysis::updatespike(st_in* in, st_out* out)
{

	//! The normal update of all neurons at spike reception.

	net->evolve_spike(calcJacobian);

	out->SC++;

	if (in->train.size() > 0)
		calc_spikeTrain(in, out);

	if (in->LyapunovExponents > 0)
		calc_LyapunovExponents(out);

	if (in->ISIstats > 0)
		calc_ISIstatistics(in, out);

	//! Save the phases of the neurons at the provided spikes in the network.
	if (in->phaseSpikes.size())
		if ((phaseSpikeCount < in->phaseSpikes.size()) && (out->SC == in->phaseSpikes[phaseSpikeCount]))
		{
			out->phaseNeurons[phaseSpikeCount] = net->get_phase_Nloc();
			if (myID == 0) out->phaseTimes[phaseSpikeCount] = out->TC;

			//! go to the next provided spike number
			phaseSpikeCount++;
		}
	
	//! Add the spike to the appropriate bin
	if ( in->instPopRateBinSize > 0 )
		out->instPopRate[int(out->TC / in->instPopRateBinSize)] += 1;
}



void LE_analysis::postprocessing(st_in* in, st_out* out)
{
	if (myID == 0)
	{
		cout << endl << "postprocessing ... " << endl;

		if (in->LyapunovExponents)
		{
			post_LyapunovExponents(in, out);

			if (in->CLV)
				calc_LyapunovVectors(out);
		}

		if (in->ISIstats > 0)
			post_ISIstatistics(in, out);

		//! normalize the instantaneous population spike rate by bin szie and neuron number
		if (in->instPopRateBinSize > 0)
			for(unsigned n=0; n<out->instPopRate.size(); n++)
				out->instPopRate[n] /= in->instPopRateBinSize*net->get_N();
	}

	//! Save the final state of the network.
	if (in->saveFinalState)
	{
		cout << "saving the final state and external currents of the network ..." << endl;

		//! Get the currents of all N neurons.
		out->finalCurrents = net->get_externalCurrents_Nloc();

		//! Get the states of the neurons.
		out->finalStates = net->get_state_Nloc();
	}
	
}



void LE_analysis::setRate_warmup(st_in* in, st_out* out)
{
	//! setRate_warmup finds the wanted rate and does a warmup of the network.

	//! if wanted and possible, find the right external currents to yield a specific average firing rate\a rateWnt
	if (in->rateWnt > 0)
	{
		if (myID == 0) cout << "adapt the external currents to yield the wanted average firing rate"<< endl;
		if (myID == 0) cout << "rateWnt = " << in->rateWnt*1000 << " Hz with precision pR = " << in->pR << " and SR = " << in->SR << " or TR = " << in->TR << " ms ..." << endl;

		bisection_externalCurrent(in, out);
	}
	else
		if (myID == 0) cout << "The provided external currents Iext will be used for the simulation." << endl;

	//! run the warmup
	if (myID == 0) cout << endl << "warmup with SW = " << in->SW << " or TW = " << in->TW << " ms ..." << endl;

	out->TW = in->TW;
	out->SW = in->SW;
	simple_iterations(&out->SW, &out->TW);

	out->rateW = out->SW/out->TW/net->get_N();
	if (myID == 0) cout << "\tspikes: " << out->SW << "\ttime: " << out->TW/1000 << "s\t -> avg. rate: " << out->rateW*1000 << " Hz" << endl;

}

void LE_analysis::simple_iterations(long long* spikes, reell* time)
{

	//! simple_iterations simply iterates the network from the current state for at least\a spikes spikes and \atime time.
	//! The given arguments\aspikes and\a time will be overwritten by the actual number of spikes and time of these iterations.

	long long spikesMin = *spikes;
	reell timeMin = *time;

	*spikes = 0;
	*time = 0;

	while (!((*spikes >= spikesMin) && (*time >= timeMin)))
	{
		//! run one iteration of the network
		reell dt = net->find_nextSpikes();
		net->evolve_dt(dt);
		net->evolve_spike(false);

		*time += dt;
		*spikes += 1;
	}
}

void LE_analysis::bisection_externalCurrent(st_in* in, st_out* out)
{
	//! find the external currents to match the wanted firing rate by bisection, the maximal number of iterations is 30
	int iterMax = 30;			//iterMax is OK to be hard coded I think
	vector<double> Itmp = in->Iext;

	reell rateTmp = 0, rateUp = 0, rateDown = 0;
	reell factorTmp = 1, factorUp = 1, factorDown = 0;

	//! run the bisection with linear guess, since the firing rate is linear to the external current in the balanced state
	bool calc = true;
	int iter = 0;

	while ((calc) && (iter < iterMax))
	{
		if (in->homogNet)
			in->Iext[0] = factorTmp*Itmp[0];
		else
			for (int n=0; n<net->get_Nloc(); n++)
				in->Iext[n] = factorTmp*Itmp[n];

		net->set_externalCurrents_Nloc(in);
		net->set_state_Nloc(in->init);				//!The reset to the initial conditions might not be too important here and could be dropped for a speed up in large networks.

		long long spikesTmp = in->SR;
		reell timeTmp = in->TR;
		simple_iterations(&spikesTmp, &timeTmp);

		rateTmp = spikesTmp/timeTmp/net->get_N();

		if (myID == 0) cout << "\t" << factorTmp << "*Iext yielded f = " << rateTmp*1000 << " Hz" << endl;

		calc = (abs(1 - rateTmp/in->rateWnt) > in->pR);

		if (calc)
		{
			if (rateTmp > in->rateWnt)
			{
				factorUp = factorTmp;
				rateUp = rateTmp;
			}
			else
			{
				factorDown = factorTmp;
				rateDown = rateTmp;
			}

			//! guess the next factor for the next current (assuming a linear slope)
			factorTmp = factorDown + (in->rateWnt - rateDown)/(rateUp - rateDown)*(factorUp - factorDown);

			//! At the beginning the current rate can be below the wanted rate, then the following is necessary to find the initial upper bound for the bisection.
			if (rateUp < in->rateWnt)
			{
				factorUp *= 2;
				factorTmp = factorUp;
			}

		}
		iter++;
	}

	if (calc)
		if (myID == 0)
			{
			cout << "the external current couldn't be found to yield the desired firing rate rateWnt = ";
			cout << in->rateWnt*1000 << " Hz with precision pR = " << in->pR << endl;
			}

	if (myID == 0) cout << "the external currents are set to " << factorTmp << "*Iext" << endl;
}



void LE_analysis::calc_spikeTrain(st_in* in, st_out* out)
{
	if (myID == 0)
	{
		vector<st_spike> spikeTmp = net->get_spikes();

		//! Check if the spiking neurons are part of the wanted neurons for the spike train.
		//! Don't know if this could be done more efficiently than going through the train vector all the time.

		for(unsigned s=0; s<spikeTmp.size(); s++)
			for(vector<int>::iterator t=in->train.begin(); t<in->train.end(); t++)
				if (spikeTmp[s].neuron == *t)
				{
					//!replace the time in spikeTmp, which was length of this interval, with the actual time\a time
					spikeTmp[s].time = out->TC;

					//!< add the currently spiking neurons\a spikeTmp to the output spike vector\a out.spikes
					out->spikeTrain.push_back(spikeTmp[s]);

					//! don't need to check for another occurence of the same neuron in\a in->train
					break;
				}
	}
}



void LE_analysis::pre_ISIstatistics(st_in* in, st_out* out)
{
	if (!(in->ISIneurons.size() > 0))
		in->ISIstats = 0;

	if (myID == 0)
	{
		//! initialize the vectors
		switch (in->ISIstats)
		{
			case 4:
				out->kurtosisNeurons = vector<reell> (in->ISIneurons.size());
				// no break, since all lower moments need to be calculated as well

			case 3:
				out->skewnessNeurons = vector<reell> (in->ISIneurons.size());
				// no break

			case 2:
				out->cvNeurons = vector<reell> (in->ISIneurons.size());
				// no break

			default:
				out->rateNeurons = vector<reell> (in->ISIneurons.size());
				spikesOut = vector<int> (in->ISIneurons.size());
				lastSpike = vector<reell> (in->ISIneurons.size());
				// no break
		}
	}
}

void LE_analysis::calc_ISIstatistics(st_in* in, st_out* out)
{
	if (myID == 0)
	{
		//! get the spike vector
		vector<st_spike> spikes = net->get_spikes();

		//! calculate the moments of the inter spike interval distribution
		for(vector<st_spike>::iterator s = spikes.begin(); s < spikes.end(); s++)
		{
			//! Check whether the spiking neuron is one of the considered ones.
			unsigned n = 0;
			while(n < in->ISIneurons.size())
			{
				if ((*s).neuron == in->ISIneurons[n])
					break;
				else
					n++;
			}

			// If the neurons wasn't found, it's not considered. Continue with the next spike.
			if (n == in->ISIneurons.size())
				continue;

			//! count the spikes to calculate the neuron's rates
			spikesOut[n]++;

			//! save the time of the first spike in the rate vector
			if (spikesOut[n] == 1)
				out->rateNeurons[n] = out->TC;

			//! calculate the isi moments
			if (in->ISIstats > 0)
			{
				if (spikesOut[n] > 0)			// if there was one spike before, we can calculate the isi
				{
					reell isi = out->TC - lastSpike[n];
					reell isi2 = isi*isi;

					switch (in->ISIstats)
					{
						case 4:
							out->kurtosisNeurons[n] += isi2*isi2;
							// no break, since all lower moments need to be calculated as well

						case 3:
							out->skewnessNeurons[n] += isi2*isi;
							// no break
						case 2:
							out->cvNeurons[n] += isi2;
							// no break
					}
				}

				lastSpike[n] = out->TC;
			}
		}
	}
}

void LE_analysis::post_ISIstatistics(st_in* in, st_out* out)
{

	for (unsigned n=0; n<in->ISIneurons.size(); n++)
		//! Calculate the statistics only if there were more than 6 spikes.
		if (spikesOut[n] > 6)
		{
			//! Calculate the rates.

			// out->rateNeurons[n] currently holds the time of the first spike
			// lastSpike[n] currently holds the time of the last spike
			out->rateNeurons[n] = (spikesOut[n] - 1)/(lastSpike[n] - out->rateNeurons[n]);


			// Calculate the moments by dividing with the number of spikes of this neuron which where considered (the first spike wasn't considered => -1).
			switch (in->ISIstats)
			{
				case 4 :
					out->kurtosisNeurons[n] /= spikesOut[n] - 1;
					// no break, since all lower moments need to be calculated as well

				case 3 :
					out->skewnessNeurons[n] /= spikesOut[n] - 1;
					//no break

				case 2 :
					out->cvNeurons[n] /= spikesOut[n] - 1;
					//no break
			}

			switch (in->ISIstats)
			{
				case 4 :
					//! Calculate the 4th standardized moment (kurtosis)
					out->kurtosisNeurons[n] = out->kurtosisNeurons[n] -
											4*out->skewnessNeurons[n]/out->rateNeurons[n] +
											6*out->cvNeurons[n]/gsl_pow_2(out->rateNeurons[n]) -
											3/gsl_pow_4(out->rateNeurons[n]);

					// standardize by divding with the 4th power of the standard deviation
					out->kurtosisNeurons[n] /= gsl_pow_2(out->cvNeurons[n] - gsl_pow_2(1./out->rateNeurons[n]));
					// no break, since all lower moments need to be calculated as well

				case 3 :
					//! Calculate the 3rd standardized moment (skewness)
					out->skewnessNeurons[n] += -3*out->cvNeurons[n]/out->rateNeurons[n] + 2/gsl_pow_3(out->rateNeurons[n]);

					// standardize by divding with the 3rd power of the standard deviation
					out->skewnessNeurons[n] /= pow(out->cvNeurons[n] - gsl_pow_2(1./out->rateNeurons[n]), (reell)3./2);
					//no break

				case 2 :
					//! Calculate the coefficient of variation
					out->cvNeurons[n] = sqrt(out->cvNeurons[n]*gsl_pow_2(out->rateNeurons[n]) - 1);
					//no break
			}

		}
		else
		{
			switch (in->ISIstats)
			{
				case 4 :
					out->kurtosisNeurons[n] = GSL_NAN;
					// no break, since all lower moments need to be calculated as well

				case 3:
					out->skewnessNeurons[n] = GSL_NAN;
					// no break

				case 2 :
					out->cvNeurons[n] = GSL_NAN;
					// no break

				default :
					out->rateNeurons[n] = GSL_NAN;
					// no break
			}
		}


	//calculate the distributions
	if (in->ISIbins > 0)
		switch (in->ISIstats)
		{
			case 4 :
				out->kurtosisDist = histogram_uniform(out->kurtosisNeurons, in->ISIbins);
				// no break, since all lower moments need to be calculated as well

			case 3 :
				out->skewnessDist = histogram_uniform(out->skewnessNeurons, in->ISIbins);
				// no break

			case 2 :
				out->cvDist = histogram_uniform(out->cvNeurons, in->ISIbins);
				// no break

			default :
				out->rateDist = histogram_uniform(out->rateNeurons, in->ISIbins);
				// no break
		}

}



void LE_analysis::pre_perturbation(st_in* in, st_out* out)
{
	if (in->pertSize)
	{
		//! evolve all neurons until after a spike
		net->evolve_dt(net->find_nextSpikes());
		net->evolve_spike(false);

		if (!applyPerturbation)
		{
			if (myID == 0) cout << endl << "saving all neurons states ... " << endl;

			//! save the original state
			stateOriginal = net->get_state_Nloc();
			IextOriginal = net->get_externalCurrents_Nloc();
		}
		else
		{
			if (myID == 0) cout << endl << "applying the perturbation of size " << in->pertSize << " ... " << endl;

			//! reset to the original state
			net->set_state_Nloc(stateOriginal);
			net->set_externalCurrents_Nloc(IextOriginal);


			//! Apply the phase perturbations
			//! This only works with the phase models so far.
			vector<reell> phasePerturbed = net->get_phase_Nloc();;

			for(unsigned n=0; n<phasePerturbed.size(); n++)
				phasePerturbed[n] += in->pertSize*in->pertVector[n];

			//!set the perturbed phase
			//! if a perturbation sets the phase above the threshold, the phase is set to the threshold in set_phase
			net->set_phase_Nloc(phasePerturbed);
		}
	}

	if (in->pertSpike || in->pertSynapse)
	{
		//! evolve all neurons until the next spike
		net->evolve_dt(net->find_nextSpikes());

		if (!applyPerturbation)
		{
			if (myID == 0) cout << endl << "saving all neurons states ... " << endl;

			//! save the original state after the spike transmission
			stateOriginal = net->get_state_Nloc();
			IextOriginal = net->get_externalCurrents_Nloc();

			//! Let the next neuron spike
			net->evolve_spike(false);

		}
		else
		{
			//!reset to the original state
			net->set_state_Nloc(stateOriginal);
			net->set_externalCurrents_Nloc(IextOriginal);


			//! Don't let the neuron spike but reset the spiking neuron.

			//! Get the next spiking neuron and save its synapses.
			vector<st_spike> spikes = net->get_spikes();
			vector<st_synapse> postOrg = in->synapses[spikes[0].neuron];

			if (in->pertSpike)
			{
				if (myID == 0) cout << endl << "skipping one spike ... " << endl;

				//! Delete all synapses of this neuron temporarily.
				in->synapses[spikes[0].neuron] = vector<st_synapse> (0);
			}
			else //if (in->pertSynapse)
			{
				if (myID == 0) cout << endl << "skipping one synaptic transmission ... " << endl;

				//! Delete the first synapse of this neuron temporarily.
				in->synapses[spikes[0].neuron].erase(in->synapses[spikes[0].neuron].begin());
			}

			//! Run the spike update. Since there are no postsynaptic neurons for the spiking neuron stored, this spike will be skipped, but the spiking neuron reset.
			net->evolve_spike(false);

			//! Reset the postsynaptic connections to the original state.
			in->synapses[spikes[0].neuron] = postOrg;
		}
	}

	//! Prepare stuff to measure phases
	phaseTimeCount = 0;
	phaseSpikeCount = 0;

	if (in->phaseSpikes.size())
	{
		in->phaseTimes.clear();
		out->phaseNeurons.resize(in->phaseSpikes.size());
		if (myID == 0) out->phaseTimes.resize(in->phaseSpikes.size());	//times need only be saved on root
	}
	else if (in->phaseTimes.size())
	{
		out->phaseNeurons.resize(in->phaseTimes.size());
		if (myID == 0) out->phaseTimes.resize(in->phaseTimes.size());	//times need only be saved on root
	}
	else
	{
		out->phaseNeurons.clear();
		out->phaseTimes.clear();
	}

	//! Save the first phases if phaseSpikes == 0.
	if ((in->phaseSpikes.size()) && (in->phaseSpikes[phaseSpikeCount] == 0))
	{
		out->phaseNeurons[phaseSpikeCount] = net->get_phase_Nloc();
		if (myID == 0) out->phaseTimes[phaseSpikeCount] = out->TC;

		//! go to the next provided spike number
		if (phaseSpikeCount < in->phaseSpikes.size() - 1)
			phaseSpikeCount++;
	}
}



void LE_analysis::warmup_ONS(st_out* out)
{
	//! Find a good step size for the orthonormalizations and warmup the ONS with SW spikes


	reell condNo = 0;
	reell condMax = 42;			// 8-)

	int NLE = out->LyapunovExponentsONS.size();

	if ((NLE > 1) && (out->ONstep > 1))
	{
		if (myID ==0 )
			cout << endl << "optimize ON step size (maximal conditional number = " << condMax << ")..." << endl;

		do
		{
			//! Find the local ONS by reiterating the orthonormalization three times
			reell dt = net->find_nextSpikes();
			net->evolve_dt(dt);
			net->evolve_spike(true);

			for(int iter=0; iter<3; iter++)
			{
				net->multiply_JacobianONS(ons, NLE);
				ons->orthonormalizeONS(NLE);
			}

			//! Evaluate the condition number of the Jacobian matrix from 6 calculations with a certain step size
			//! If the condition number is too large (here larger than 42) than the steps size needs to be smaller
			condNo = 0;
			for(int bla=0; bla<6; bla++)
			{
				for(int st=0; st<out->ONstep; st++)
				{
					reell dt = net->find_nextSpikes();
					net->evolve_dt(dt);
					net->evolve_spike(true);

					net->multiply_JacobianONS(ons, NLE);
				}

				ons->orthonormalizeONS(NLE);

				condNo += ons->get_normONS(0)/ons->get_normONS(NLE-1);
			}
			condNo /= 6;

			if (myID == 0) cout << "\tthe average condition number with ON steps = " << out->ONstep << " is " << condNo << endl;

			if (condNo > condMax)
				out->ONstep = (int)floor(out->ONstep/2);

		}
		while ((out->ONstep > 1) && (condNo > condMax));

		if (myID ==0 ) cout << "\tsetting ON step size to " << out->ONstep << endl;
	}



	//! warmup of the ONS with the found step size with 1 spike per neuron
	if (myID == 0) cout << endl << "warmup of the orthonormal system with " << out->SWONS << " spikes ..." << endl;

	int long long s = 0;
	while (s < out->SWONS)
	{
		for (int st=0; st<out->ONstep; st++)
		{
			reell dt = net->find_nextSpikes();
			net->evolve_dt(dt);
			net->evolve_spike(true);

			net->multiply_JacobianONS(ons, NLE);

			s++;
		}

		ons->orthonormalizeONS(NLE);
	}

}

void LE_analysis::pre_LyapunovExponents(st_in* in, st_out* out)
{
	//! initialize the orthonormal system
	//! The state dimension is the same for all neurons, since net->get_allNeuronSame must be true for this calculation. This is checked above.
	delete ons;
	ons = new LE_ons(in->LyapunovExponents, net->get_Nloc(), net->get_stateDim(0), in->seedONS + myID);

	//! zero vector for Lyapunov exponents in the out structure (on all nodes including the slaves, because NLE is extract from LyapunovExponents.size() later on)
	out->LyapunovExponentsONS = vector<reell> (in->LyapunovExponents);

	//! initial values for the convergence calculation
	LEmaxOld = 0;
	LEminOld = 0;
	out->pLEONS = in->pLEONS;

	//! the orthonormalization step size will be adapted and the ONS warmuped
	out->ONstep = (in->ONstep > 0) ? in->ONstep : 1;
	out->SWONS = in->SWONS;

	if (in->LyapunovExponents > 0)
		warmup_ONS(out);

	//! Important: The result during the warmup are excluded in the following calculations including the convergence vector.

	//! Should the convergence of the Lyapunov exponents be saved (necessary only on root)?
	if ((myID == 0) && in->LyapunovExponentsConvergence)
	{
		//! initialize and assign the first output vector element with 0
		out->LEconvergence = vector<vector<reell> > (1, out->LyapunovExponentsONS);

		//! Reserve the minimal amount of memory for the number of orthonormalizations in case the number of spikes per neuron were provided.
		//! If the time was provided this can't be done due to the lack of knowledge of the minimal number of orthonormalizations.
		if (in->SC)
			out->LEconvergence.reserve(in->SC/out->ONstep);

	}
	else
	{
		out->LEconvergence = vector<vector<reell> > (0);
		out->LEtimes = vector<reell> (0);
	}

	//! Store the times of orthonormalization if the convergence is saved or the covariant Lyapunov exponents are computed
	if ((myID == 0) && (in->LyapunovExponentsConvergence || in->CLV))
	{
		out->LEtimes = vector<reell> (1);
		if (in->SC)
			out->LEtimes.reserve(in->SC/out->ONstep);
	}


	//! If the covariant Lyapunov vectors are to be calculated the projection matrix R needs to be stored for each orthonormalization.
	//! This should happen after the warmup of the ons, since this should not be used for the covariant Lyapunov vector calculation.
	if ((myID == 0) && in->CLV)
	{
		ons->saveProjections = true;

		//! Initialize the covariant Lyapunov vectors.
		ons->initCLV();

		out->SWCLV = in->SWCLV;
		out->LyapunovExponentsCLV = vector<reell> (in->LyapunovExponents);

		out->localLyapunovExponents = vector<vector<reell> > (1, out->LyapunovExponentsCLV);
		if (in->SC)
			out->localLyapunovExponents.reserve((in->SC - in->SWCLV)/out->ONstep);
	}
	else
		ons->saveProjections = false;

}

void LE_analysis::calc_LyapunovExponents(st_out* out)
{
	int NLE = out->LyapunovExponentsONS.size();
	net->multiply_JacobianONS(ons, NLE);

	//! Reorthonormalize only every other step to save computation time
	if(out->SC % out->ONstep == 0)
	{
		ons->orthonormalizeONS(NLE);

		for(int n=0; n<NLE; n++)
			if (myID == 0) out->LyapunovExponentsONS[n] += log(ons->get_normONS(n));


		//! Should the convergence be monitored?
		if ((myID == 0) && (out->LEtimes.size()))
		{
			//save the current times
			out->LEtimes.push_back(out->TC);

			//save the accumulated log(norms) (divide by the time in postprocessing, then the element index is save to use)
			out->LEconvergence.push_back(out->LyapunovExponentsONS);
		}


		//! Check the convergence of the largest and the smallest Lyapunov exponent (3 ONsteps back in time)
		if ((out->pLEONS > 0) && (int(out->SC/out->ONstep) % 3 == 0))
		{
			reell LEmax = out->LyapunovExponentsONS[0]/out->TC;
			reell LEmin = out->LyapunovExponentsONS[NLE-1]/out->TC;

			bool notConvMax = (fabs(LEmax) > 1e-10) ? (fabs(1 - LEmaxOld/LEmax) > out->pLEONS) : false;
			bool notConvMin = (fabs(LEmin) > 1e-10) ? (fabs(1 - LEminOld/LEmin) > out->pLEONS) : false;

			notConverged = notConvMax || notConvMin;

			LEmaxOld = LEmax;
			LEminOld = LEmin;
		}
	}
}

void LE_analysis::post_LyapunovExponents(st_in* in, st_out* out)
{
	for(unsigned n=0; n<out->LyapunovExponentsONS.size(); n++)
		out->LyapunovExponentsONS[n] /= out->TC;

	for (unsigned t=1; t<out->LEconvergence.size(); t++)
	{
		// divide the stored log(norms) by the time to yield the Lyapunov exponents (not the first element, it's zero)
		for(unsigned n=0; n<out->LEconvergence[t].size(); n++)
			out->LEconvergence[t][n] /= out->LEtimes[t];
	}
}

void LE_analysis::calc_LyapunovVectors(st_out* out)
{
	if (myID == 0)
	{
		//! Get the number of orthonormalizations for which the projections were stored.
		unsigned SC = ons->get_numberSavedProjections();
		int NLE = out->LyapunovExponentsONS.size();
		unsigned SW = out->SWCLV/out->ONstep;

		cout << "calculating the covariant Lyapunov vectors from " << SC << " orthonormalization steps (this always runs on only 1 processor) ..." << endl;

		int Sdone = 1;
		//! Run the covariant Lyapunov vector calculation.
		for (unsigned s=1; s<=SC; s++)
		{

			ons->backwardIterationCLV(SC-s, NLE);

			//! After the warmup of the CLVs start the analysis
			if (s > SW)
			{
				//! Get the log norms of the current orthonormalization of the CLVs
				vector<reell> logNorms(NLE);
				for(int n=0; n<NLE; n++)
					logNorms[n] = -log(ons->get_normCLV(n));		// -(minus) because time is reversed

				//! store the local Lyapunov exponents.
				//! They are in reverse order, the 0 element is the first after the warmup of the CLVs. Change order later.
				out->localLyapunovExponents.push_back(logNorms);

				//! store the backward Lyapunov exponents.
				for(int n=0; n<NLE; n++)
					out->LyapunovExponentsCLV[n] += logNorms[n];
			}

			if (s == SW)
				cout << "\twarmup done, step: " << s << endl;

			if  (s >= SC/10.*Sdone)
				if (myID == 0) cout << "\t" << Sdone++*10 << "% done ... step: " << s << endl;


		}

		//! Postprocessing of the results

		//! Calculate the backward LEs.
		reell timeCLV = out->LEtimes[SC - (SW+1)] - out->LEtimes[0];
		for(int n=0; n<NLE; n++)
			out->LyapunovExponentsCLV[n] /= timeCLV;

		//! Swap order of the local LEs.
		for (unsigned s=0; s<out->localLyapunovExponents.size()/2; s++)
			out->localLyapunovExponents[s].swap(out->localLyapunovExponents[out->localLyapunovExponents.size() - (s+1)]);

	}
	else
	{
		cout << "The covariant Lyapunov vector calculation is completely done on the root process, but was called on a slave here" << endl;
		throw(2);
	}
}



vector<vector<reell> > LE_analysis::histogram_uniform(vector<reell> data, int bins)
{
	//! Calculates the histogram with uniformly spaced bins on the x axis, the x-values return the center of the bins.

	//! Erase NANs from the distribution
	unsigned x = 0;
	do
	{
		if (gsl_isnan(data[x]))
			data.erase(data.begin() + x);
		else
			x++;
	}
	while (x < data.size());

	//! sort data in ascending order
	sort(data.begin(), data.end());

	//! set the uniformly distributed x values
	vector<vector<reell> > hist(bins, vector<reell> (2));

	reell dataMin = data[0];
	reell dataMax = data[data.size() - 1];

	reell dx = (dataMax - dataMin)/bins;
	reell dx2 = dx/2;

	for(int b=0; b<bins; b++)
	{
		hist[b][0] = dataMin + b*dx + dx2;					//centers of the bins
		hist[b][1] = 0;
	}

	//! fill the histogram
	vector<reell>::iterator it = data.begin();

	for(int b=0; b<bins; b++)
		while (*it <= hist[b][0] + dx2)
		{
			hist[b][1] = hist[b][1] + 1;

			if (it < data.end() - 1)
				it++;
			else
				break;			//shouldn't occur, since hist[b][bins-1] == maxData
		}

	//! normalize the histogram
	for(int b=0; b<bins; b++)
		hist[b][1] /= dx*data.size();				//all bins are equally space here (dx)!

	return hist;
}

/** possible extensions:
 *  logarithmic histograms
 */
