/*
 * LENetwork.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#include "LEnetwork.h"

LE_network::LE_network(st_in* in)
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

	gsl_rng_env_setup();
	const gsl_rng_type *TYPE = gsl_rng_mt19937;
	rngSyn = gsl_rng_alloc(TYPE);
	set_rngSynapses(1);


	//! Number of neurons in the network.
	Nall = in->Nall;
	Nloc = in->Nloc;
	Noffset = in->Noffset;




	neurons = vector<LE_neuron*> (Nloc);
	for (int n=0; n<Nloc; n++)
	{
		//! initialize the neurons
		//! if homogeneous network, then initialize all neurons like neuron 0
		int nHomo = (in->homogNet) ? 0 : n;

		if (in->neuronType[nHomo] == 1)
			neurons[n] = new LE_rapidtheta(in->rapidness[nHomo], in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 2)
			neurons[n] = new LE_lif(in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 3)
			neurons[n] = new LE_type1type2(in->type1type2Para[nHomo], in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 4)
			neurons[n] = new LE_rapidtheta(in->rapidness[nHomo], in->reset[nHomo], in->threshold[nHomo], in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 10)
			neurons[n] = new LE_twoDlinear(in->tauM[nHomo], in->twoDlinearParas[nHomo]);

		if (in->neuronType[nHomo] == 11)
			neurons[n] = new LE_clif_fast3(in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 12)
			neurons[n] = new LE_clif_fast2(in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 13)
			neurons[n] = new LE_clif_slow2(in->tauM[nHomo]);

		if (in->neuronType[nHomo] == 14)
			neurons[n] = new LE_clif_slow3(in->tauM[nHomo]);
	}

	set_externalCurrents_Nloc(in);
	set_state_Nloc(in->init);
	set_synapses(in);

	//! Test whether all neurons are of the same kind, e.g. phase neurons.
	LE_neuron *test;

	//! Test for phase neurons. They can be different phase neurons though.
	allPhaseNeurons = true;
	allClifNeurons = true;

	for (int n=0; n<Nloc; n++)
	{
		test = dynamic_cast<LE_phaseneuron*>(neurons[n]);
		if (test == 0)					//this one isn't a phase neuron
		{
			//can't calcualte the Jacobian of this network
			allPhaseNeurons = false;
			break;
		}

		test = dynamic_cast<LE_clif*>(neurons[n]);
		if (test == 0)					//this one isn't a cLIF neuron
		{
			//can't calcualte the Jacobian of this network
			allClifNeurons = false;
			break;
		}

	}

#ifdef PAR
	bool global;

	MPI_Allreduce(&allPhaseNeurons, &global, 1, MPI_SHORT, MPI_BAND, MPI_COMM_WORLD);
	allPhaseNeurons = global;

	MPI_Allreduce(&allClifNeurons, &global, 1, MPI_SHORT, MPI_BAND, MPI_COMM_WORLD);
	allClifNeurons = global;

#endif

	allNeuronsSame = allPhaseNeurons || allClifNeurons; // add other neuron types
	if (allNeuronsSame == false)
		cout << "-------------" << "Not all the neurons are phase neurons. Thus, the single spike Jacobian cannot be calculated." << "-------------" << endl;

	//! Initilize the Jacobian vector construct.
	//! the 4th dimension -> This will have the dimension of the state variable of the neuron class squared.
	vector<reell> stateTmp(neurons[0]->get_stateDim()*neurons[0]->get_stateDim());
	//! the 3rd dimension -> This will have two dimensions, one for the diagonal and one for off-diagonal Jacobian term.
	vector<vector<reell> > jacTmp(2, stateTmp);
	//! the 2nd dimension -> This will be the number of neurons in the network (for now, set it 1, to save memory in case the Jacobian isn't calculated)
	vector<vector<vector<reell> > > neurTmp(1, jacTmp);
	//! the 1st dimension -> This will have the dimension of the number of synchronous spikes in one iteratino (for now, set it to 1).
	vector<vector<vector<vector<reell> > > > spkTmp(1, neurTmp);
	//! the Jacobian vector construct is resized to the appropriate dimensions in calc_Jacobian().
	jacobian = spkTmp;


}

LE_network::~LE_network()
{
	gsl_rng_free(rngSyn);																	//! free the random number generator

	for (int n=0; n<Nloc; n++)
		delete neurons[n];
}


//************
// set stuff *
//************


void LE_network::set_rngSynapses(unsigned seed)
{
	gsl_rng_set(rngSyn, 1 + myID + seed);
}

void LE_network::set_externalCurrents_Nloc(st_in* in)
{
	//! set the external current
	for (int n=0; n<Nloc; n++)
	{
		int nHomo = (in->homogNet) ? 0 : n;
		neurons[n]->set_externalCurrent(in->Iext[nHomo]);
	}
}

void LE_network::set_externalCurrents_Nloc(vector<reell>& currents)
{
	//! set the external current
	for (int n=0; n<Nloc; n++)
		neurons[n]->set_externalCurrent(currents[n]);

}

void LE_network::set_synapses(st_in* in)
{
	//! Set the postsynaptic matrix, we'll use a pointer to the in structure.
	//! This means that the synapse could in principle be changed dynamically within the program by just changing the in structure accordingly.
	synapses = &in->synapses;

	//! Delete autapses in the adjacency matrix, if they occur.
	//! If autapses are really wnated then the update/reset function and the single spike Jacobian need to be changed correctly.

	for (unsigned s=0; s<synapses->size(); s++)
	{
		int localID = global2local(s);

		for(vector<st_synapse>::iterator p=(*synapses)[s].begin(); p<(*synapses)[s].end(); )
			if ((*p).post == localID)
			{
				if (myID == 0) cout << "deleting autapse on neuron " << s << endl;
				(*synapses)[s].erase(p);

				//don't increase the iterator p, because one element has been removed
			}
			else
				p++;	//increase the iterator p, to check the next element
	}


/*
	for (unsigned s=0; s<synapses.size(); s++)
	{
		cout << myID << " : " << s << " = " << synapses[s].size() << " : ";
		for (unsigned m=0; m<synapses[s].size(); m++)
			cout << synapses[s][m].post << "(J=" << synapses[s][m].cplg << ", p=" << synapses[s][m].prob<< ") ";
		cout << endl;
	}
*/
}

void LE_network::set_state_Nloc(vector<vector<reell> >& states)
{
	for (int n=0; n<Nloc; n++)
		neurons[n]->set_state(states[n]);
}

void LE_network::set_phase_Nloc(vector<reell>& phase)
{
	for (int n=0; n<Nloc; n++)
		neurons[n]->set_phase(phase[n]);
}

void LE_network::set_externalCurrents_neurons(vector<int>& curNeurons, vector<reell>& currents)
{
	for (unsigned n=0; n<curNeurons.size(); n++)
		neurons[curNeurons[n]]->set_externalCurrent(currents[n]);
}


//************
// get stuff *
//************


vector<st_spike> LE_network::get_spikes()
{
	return spikes;
}

int LE_network::get_N()
{
	return Nall;
}

int LE_network::get_Nloc()
{
	return Nloc;
}

bool LE_network::get_allNeuronSame()
{
	return allNeuronsSame;
}

int LE_network::get_stateDim(int n)
{
	return neurons[n]->get_stateDim();
}

vector<vector<reell> > LE_network::get_state_Nloc()
{
	//!Returns all neurons states.

	vector<vector<reell> > state(Nloc);

	for (int n=0; n<Nloc; n++)
		state[n] = neurons[n]->get_state();

	return state;
}

vector<reell> LE_network::get_phase_Nloc()
{
	//!Returns all neurons phases. This must be correctly implemented for the 2D models (cLIF and twoDlinear)!

	vector<reell> phases(Nloc);

	for (int n=0; n<Nloc; n++)
		phases[n] = neurons[n]->get_phase();

	return phases;
}

vector<reell> LE_network::get_externalCurrents_Nloc()
{
	//! get the external current
	vector<reell> currents(Nloc);
	for (int n=0; n<Nloc; n++)
		currents[n] = neurons[n]->get_externalCurrent();

	return currents;
}

vector<reell> LE_network::get_externalCurrents_neurons(vector<int>& curNeurons)
{
	//! returns the external current of the neurons given in vector\a curNeurons
	vector<reell> current(curNeurons.size());

	for (unsigned n=0; n<curNeurons.size(); n++)
		current[n] = neurons[curNeurons[n]]->get_externalCurrent();

	return current;
}


//******************
// transform stuff *
//******************


//! transform from global neuron ID to the local neuron ID
inline int LE_network::global2local(int globalID)
{
#ifdef PAR
	return globalID - Noffset;
#else
	return globalID;
#endif
}

//! transform from global neuron ID to the local neuron ID
inline int LE_network::local2global(int localID)
{
#ifdef PAR
	return localID + Noffset;
#else
	return localID;
#endif
}

inline int LE_network::spikingNode(int globalID)
{
	//! calculate the processor to which the spiking neurons belongs in parallel mode (integer division)
#ifdef PAR

	int spikingID = globalID/(Nall/nP);

	//! If the number of all neurons cannot be divided by the number of processors, Nloc is uneven and the last
	//! node processes more neuron than the other nodes.
	if (spikingID >= nP)
		spikingID = nP - 1;

	return spikingID;

#else

	return 0;

#endif

}


//******************
// iteration stuff *
//******************


reell LE_network::find_nextSpikes()
{
	//! Calculate the next spikes and stores them in spikes, multiple ones if there are synchronous spikes.


	//! find the next spike time locally (of the neurons on each processor)
	st_spike spikeTmp;

	reell dtTmp = neurons[0]->get_spikeTime();

	spikeTmp.time = dtTmp;
	spikeTmp.neuron = 0;													//! store the global ID of the neuron

	spikes = vector<st_spike> (1, spikeTmp);											//!< fill temporarily with the first neuron's values

	//!find the next spike time
	for (int n=1; n<Nloc; n++)
	{
		reell dt = neurons[n]->get_spikeTime();											//!< get the neuron's spike time

		if (dt <= dtTmp)																//!< this neurons spikes before or with the previously selected ones
		{
			spikeTmp.time = dt;
			spikeTmp.neuron = n;

			if (dt < dtTmp)
			{
				spikes = vector<st_spike> (1, spikeTmp);								//!< overwrite the spikes vector because if this time is earlier than all previous ones.
				dtTmp = dt;
			}
			else
				spikes.push_back(spikeTmp);												//!< append this time to the synchronous spike vector.
		}
	}

#ifdef PAR

	//! 0) change to global index of the psiking neurons
	for (unsigned s=0; s<spikes.size(); s++)
		spikes[s].neuron = local2global(spikes[s].neuron);



	//! 1) Gather all spike times and the number of synchronous spikes from all processors.

	st_double_int local;
	local.d = spikes[0].time;
	local.i = spikes.size();

	vector<st_double_int> global(nP);

	MPI_Allgather(&local, 1, MPI_DOUBLE_INT, &global.front(), 1, MPI_DOUBLE_INT, MPI_COMM_WORLD);

	//! 2) Find the shortest next spike time.

	reell dt = global[0].d;

	for (int p=1; p<nP; p++)
	{
		if (global[p].d < dt)
			dt = global[p].d;
	}

	//! 3) Gather all synchronous spikes (possibly from different processors)

	vector<st_spike> spikesLocal, spikesGlobal;

	for (int p=0; p<nP; p++)
		if (global[p].d == dt)
		{
			if (myID == p)
				spikesLocal = spikes;
			else
				spikesLocal = vector<st_spike>(global[p].i);		//!< allocate the memory for broadcast

			MPI_Bcast(&spikesLocal.front(), global[p].i, MPI_DOUBLE_INT, p, MPI_COMM_WORLD);

			spikesGlobal.insert(spikesGlobal.end(), spikesLocal.begin(), spikesLocal.end());
		}

	spikes = spikesGlobal;

#endif

	return spikes[0].time;

}

void LE_network::evolve_dt(reell dt)
{
	//!This evolves all neurons states for the time interval \adt, e.g. just before the next spike.
	if (dt > 0)
	{
		//! Calculate some dummy values if necessary on neuron with index 0
		neurons[0]->dummy_evolve(dt, &dummyEvolve);

		//! Evolve all neurons.
		for (int n=0; n<Nloc; n++)
			neurons[n]->evolve_dt(dt, &dummyEvolve);
	}
}

void LE_network::evolve_spike(bool calcJacobian)
{
	//! This requires that find_nextSpikes was called before and the current spiking neurons are thus saved in \aspikes.
	//! All neurons were also evolved with evolve_dt just before the spike reception.

	/** Synchronous spikes are still treated in subsequent order according to their neuron number. This means, when two neurons
	 *  spike at the same time, the neuron with the lower neuron index will send its spike first to its postsynaptic neurons.
	 *  Then the other neuron will send its spike. In case they both share a postsynaptic neuron, the spikes arrive in subsequent
	 *  order. If the coupling strength is the same from both spiking neurons, the final state of the postsynaptic neuron is
	 *  unique, it doesn't depend on the order of the spikes to arrive. However, if the coupling strengths are different and the
	 *  phase response curve is nonlinear, then the final state of the postsynaptic neuron depends on the order of the arriving
	 *  spikes. Therefore, the network will have follow different trajectory if the indices of the neurons, e.g. that spike
	 *  synchronously, are switched without changing the network really. This could be avoided, if synchronous spikes also arrive
	 *  synchronously at their postsynaptic neurons. For the update of the neurons this would be simply result in summing the
	 *  coupling strengths, but the derivation of the single spike Jacobian would change considerably and needs to be done carefully.
	 */

	//! Calculate the single spike Jacobian (multiple in case of synchronous spikes).
	if (calcJacobian)
		//! To calculate the single spike Jacobian correctly, the previously called evolve_dt must cover the entire interval!
		calcJacobian_and_updateNeurons();
	else
		updateNeurons();
	
	//
	for (int n=0; n<Nloc; n++)
		neurons[n]->dummy_PostUpdate();
}

//! Updates all neurons that receive the current spikes without calculating the single spike Jacobian.
void LE_network::updateNeurons()
{
	//! Update the postsynaptic neurons if the release probability of the synapse permits it.
	for(vector<st_spike>::iterator s = spikes.begin(); s < spikes.end(); s++)
		for(vector<st_synapse>::iterator p = (*synapses)[(*s).neuron].begin(); p < (*synapses)[(*s).neuron].end(); p++)
			if ((*p).prob < 1)
			{
				if (gsl_rng_uniform(rngSyn) < (*p).prob)
					neurons[(*p).post]->evolve_spike((*p).cplg);
			}
			else
				neurons[(*p).post]->evolve_spike((*p).cplg);

	//! The neurons need to be reset AFTER every postsynaptic neuron was updated,
	//! in case one of the spiking neuron is one of the postsynaptic neurons of another simultaneously spiking neuron.
	for(vector<st_spike>::iterator s = spikes.begin(); s < spikes.end(); s++)
	{
#ifdef PAR

		//! check if the spiking neuron belongs to this processor for parallel mode (integer division)
		if (myID == spikingNode((*s).neuron))
			neurons[global2local((*s).neuron)]->reset();

#else

		neurons[(*s).neuron]->reset();

#endif
	}

}

//! Calculates the single spike Jacobian of phase neurons using the derivative of the phase transition curve, Eq.~(2.15)
void LE_network::calcJacobian_and_updateNeurons()
{
	/** This calculates the single spike Jacobian stored in a sparse way, since there are at most one diagonal element
	 *  and one nondiagonal element per neuron. Important is that this method is called before the update of the
	 *  postsynaptic neurons and after the evolution of all neurons for the full time interval, because the state before
	 *  the spike reception and the length of the time interval defined in each neuron class are used to calculate the
	 *  single spike Jacobian elements.
	 */

	//! Note, resize doesn't do anything if the vector is already large enough, so it's cheap, although it seems that it is done every time.
	//! Resize the 2nd dimension to the number of synchronous spikes in this interval.

	jacobian.resize(spikes.size(), jacobian[0]);

	for(unsigned s=0; s<spikes.size() ; s++)
	{
		//! Resize the 1st dimension to the number of neurons.
		jacobian[s].resize(Nloc, jacobian[s][0]);

		//!< The spiking neuron.
		int spkNeuron = spikes[s].neuron;


		//! There is one variable from the spiking neuron needed, e.g. its velocity at the spike.
#ifdef PAR
		//! calculate the processor to which the spiking neurons belongs in parallel mode
		int spkID = spikingNode(spkNeuron);
		int spkNeuronLocal = (myID == spkID) ? global2local(spkNeuron) : -1;

		if (myID == spkID)
			neurons[spkNeuronLocal]->dummy_Jacobian_spiking(&dummyJacobian);

		//! broadcast the size of the dummy variable
		int dummySz = (myID == spkID) ? dummyJacobian.size() : 0;
		MPI_Bcast(&dummySz, 1, MPI_DOUBLE, spkID, MPI_COMM_WORLD);

		//! broadcast this dummy variable to all processors
		if (myID != spkID)
			dummyJacobian.resize(dummySz);
		MPI_Bcast(&dummyJacobian.front(), dummySz, MPI_DOUBLE, spkID, MPI_COMM_WORLD);

#else
		neurons[spkNeuron]->dummy_Jacobian_spiking(&dummyJacobian);
#endif

		vector<st_synapse> noSynapses(1);			//!< dummyJacobian vector if the spiking neuron doesn't have any postsynaptic neurons
		noSynapses[0].post = -1;					//!< set this to -1 such that (n != (*p).post) in the following loop
													//!< if there are postsynaptic neurons, then everythin is fine

		//! iterator @ postsynaptic neurons
		vector<st_synapse>::iterator p = ((*synapses)[spkNeuron].size() > 0) ? (*synapses)[spkNeuron].begin() : noSynapses.begin();
		bool release = ((*p).prob < 1) ? (gsl_rng_uniform(rngSyn) < (*p).prob) : true;

		//! The Jacobian elements for all neurons.
		for (int n=0; n<Nloc; n++)
		{
			//! The 3rd dimension is the square of the state dimension.
			//! The 4th dimension contains the diagonal and nondiagonal element of the single spike Jacobian (2 elements)

			//! The Jacobian elements for the postsynaptic neurons.
			//!< The postsynaptic neurons indices in\a synapses must be in increasing order.

			if (n == (*p).post)
			{
				if (release)
				{
					jacobian[s][n] = neurons[(*p).post]->calc_JacElem_postsynaptic(&dummyJacobian, (*p).cplg);
					neurons[(*p).post]->evolve_spike((*p).cplg);
				}
				else
				{
					//! treat this neuron, as if it wouldn't be connected to the spiking neuron
					//! the diagonal elements are calculated as if this neuron was in the not-postsynaptic set
					jacobian[s][n] = neurons[n]->calc_JacElem_self(&dummyJacobian);

					//! the nondiagonal elements are set zero, but they must exist, because in the multiplication
					//! with the ONS, it is not known which synapse was deactivated here, thus for all postsynaptic
					//! neurons the diagonal and nondiagonal elements are expected and must exist
					jacobian[s][n].push_back(vector<reell> (jacobian[s][n][0].size()));
				}

				if (p < (*synapses)[spkNeuron].end() - 1)
				{
					p++;
					release = ((*p).prob < 1) ? (gsl_rng_uniform(rngSyn) < (*p).prob) : true;
				}
			}

			//! The Jacobian elements for the spiking neuron.
#ifdef PAR

			else if ((n == spkNeuronLocal) && (myID == spkID))
				jacobian[s][n] = neurons[spkNeuronLocal]->calc_JacElem_spiking(&dummyJacobian);

#else

			else if (n == spkNeuron)
				jacobian[s][n] = neurons[spkNeuron]->calc_JacElem_spiking(&dummyJacobian);

#endif

			//! The Jacobian elements for all other neurons that are not postsynaptic to the spiking neurons.
			else
				jacobian[s][n] = neurons[n]->calc_JacElem_self(&dummyJacobian);
		}
	}

	//! The neurons need to be reset AFTER every postsynaptic neuron was updated,
	//! in case one of the spiking neuron is one of the postsynaptic neurons of another simultaneously spiking neuron.
	for(vector<st_spike>::iterator s = spikes.begin(); s < spikes.end(); s++)
	{
#ifdef PAR

		//! check if the spiking neuron belongs to this processor for parallel mode
		if (myID == spikingNode((*s).neuron))
			neurons[global2local((*s).neuron)]->reset();


#else

		neurons[(*s).neuron]->reset();

#endif
	}
}

void LE_network::multiply_JacobianONS(LE_ons* ons, int NLE)
{
	if (ons->get_COL() < NLE)
	{
		cout << "The number of Lyapunov exponents " << NLE << " to be calculated is too big." << endl;
		throw(1);
	}

	if (ons->get_DIM()*ons->get_DIM() != (int)jacobian[0][0][0].size())
	{
		cout << "DIM (state) = " << ons->get_DIM() << "\tDIM (jacobian) = " << jacobian[0][0][0].size() << endl;
		cout << "The dimension of the state variable of the ons and that of the Jacobian do not match." << endl;
		cout << "The Jacobian dimension needs to be the state variable dimension squared." << endl;
		throw(1);
	}

	if (ons->get_ROW() != Nloc)
	{
		cout << "The dimension of the ons vectors " << ons->get_ROW() << " and the number of neurons " << Nloc << " do not match." << endl;
		throw(1);
	}

	int DIM  = ons->get_DIM();
	int ROW  = ons->get_ROW();

	reell onsJ[DIM];
	reell onsTmp[DIM];

	//! Mutiply the Jacobain with the ons for each synchronous spike.
	for(unsigned s=0; s<spikes.size(); s++)
	{
		//!< The spiking neuron.
		int spkNeuron = spikes[s].neuron;

#ifdef PAR

		int spkID = spikingNode(spkNeuron);
		int spkNeuronLocal = (myID == spkID) ? global2local(spkNeuron) : -1;

#endif

		for(int c=0; c<NLE; c++)
		{
			//! Save the spiking neurons elements for the calculation of the nondiagonal elements.
#ifdef PAR

			if (myID == spkID)
				for(int din=0; din<DIM; din++)
					onsJ[din] = ons->ons[c][spkNeuronLocal*DIM + din];

			//! broadcast these ons elements to all processors
			MPI_Bcast(onsJ, DIM, MPI_DOUBLE, spkID, MPI_COMM_WORLD);

#else

			for(int din=0; din<DIM; din++)
				onsJ[din] = ons->ons[c][spkNeuron*DIM + din];

#endif


			vector<st_synapse> noSynapses(1);			//!< dummy vector if the spiking neuron doesn't have any postsynaptic neurons
			noSynapses[0].post = -1;					//!< set this to -1 such that (n != (*p).post) in the following loop
														//!< if there are postsynaptic neurons, then everythin is fine

			//! iterator @ postsynaptic neurons
			vector<st_synapse>::iterator p = ((*synapses)[spkNeuron].size() > 0) ? (*synapses)[spkNeuron].begin() : noSynapses.begin();

			for(int n=0; n<ROW; n++)
			{
				for(int dout=0; dout<DIM; dout++)
					onsTmp[dout] = 0;

				//! The diagonal elements of the phase neuron models are the identity, skip this.
				for(int dout=0; dout<DIM; dout++)
					for(int din=0; din<DIM; din++)
						onsTmp[dout] += jacobian[s][n][0][dout*DIM + din]*ons->ons[c][n*DIM + din];

				//! The additional nondiagonal Jacobian elements for the postsynaptic neurons.
				//! The postsynaptic neurons indices in\a synapses must be in increasing order.
				if (n == (*p).post)
				{
					for(int dout=0; dout<DIM; dout++)
						for(int din=0; din<DIM; din++)
							onsTmp[dout] += jacobian[s][n][1][dout*DIM + din]*onsJ[din];

					if (p < (*synapses)[spikes[s].neuron].end() - 1)
						p++;
				}

				for(int dout=0; dout<DIM; dout++)
					ons->ons[c][n*DIM + dout] = onsTmp[dout];
			}
		}
	}
}


/*
inline void LE_network::init_randomgraph(reell prob, reell J)
{
	if (myID == 0) cout << "initialize random network with AVERAGE indegree K per node ..." << endl;

	//init
	(*synapses) = vector<vector<st_synapse> > (Nall);

	st_synapse synTmp;
	synTmp.post = 0;
	synTmp.cplg = J;

	for(int n=0; n<Nall; n++)
		for (int m=0; m<Nloc; m++)
			if (gsl_rng_uniform(rng) < prob)
				if (n != m)					//no autapses (self loops)
				{
					synTmp.post = m;
					(*synapses)[n].push_back(synTmp);
				}
}
*/
