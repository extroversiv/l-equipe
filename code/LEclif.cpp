/*
 * LEclif.cpp
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#include "LEclif.h"

LE_clif::LE_clif(reell tau) : LE_neuron(tau, vector<reell> (1, 0), vector<reell> (1, 1))
{
	//! The rheobase current is always 1 in the LIF model.
	//! The reset value is 0 and the threshold is 1.
	Irheo = 1;

	//! The state variable is two-dimensional in the CLIF model.
	state = vector<reell> (2);
}


LE_clif::~LE_clif()
{

}

void LE_clif::set_phase(reell&)
{
	cout << "reell LE_clif::set_phase() not implemented yet!" << endl;
	throw(2);
}

reell LE_clif::get_phase()
{
	cout << "reell LE_clif::get_phase() not implemented yet!" << endl;
	throw(2);
}

void LE_clif::set_externalCurrent(reell& Iext)
{
	this->Iext = Iext;
	calcSpikeTime = true;;
}

void LE_clif::evolve_dt(reell dt, vector<reell>* dummy)
{
	dt /= tauM;

	//! membrane time exponential and synaptic time exponential
	reell expM = exp(-dt);
	reell expS = powGamma(expM);						//!< power function depending on gamma

	reell tmpI = state[1]/(1 - gamma);

	//! evolve the state between spikes
	state[0] =  1 + Iext + tmpI*expS - (1 + Iext + tmpI - state[0])*expM;		//!< Eq. (4.8)
	state[1] *= expS;

	//! update the spike time
	calcSpikeTime = true;
	//! We cannot just use spikeTime -= dt because this acccumulates precision losses. For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.

	//! Sum up the time since the last spike in the network (actually expM and expS here to save some computations)
	jacISIexpM *= expM;
	jacISIexpS *= expS;

}

void LE_clif::evolve_spike(reell c)
{
	//! \a state[0] - the voltage doesn't change when receiving a spike, only the current\a state[1] changes
	state[1] += c*gamma;
	calcSpikeTime = true;;
}

//! Set the state variable to the reset value.
void LE_clif::reset()
{
	//! Only the voltage is reset, the current stays unaffected.
	state[0] = stateReset[0];
	calcSpikeTime = true;;
}

void LE_clif::dummy_PostUpdate()
{
	jacISIexpM = 1;
	jacISIexpS = 1;
}

//! The diagonal elements of the single spike Jacobian for the postsynaptic neurons.
vector<vector<reell> > LE_clif::calc_JacElem_postsynaptic(vector<reell>* dummy, reell c)
{
	vector<vector<reell> > jacTmp(2, vector<reell>(4));					//!< the vector with Jacobian elements
	c /= tauS;															//!< use c for gamma*c/tauM = c/tauS

	//! 4 elements for the diagonal and the nondiagonal part, each.
	jacTmp[0][0] = jacISIexpM;											//!< the diagonal element of the voltage-voltage block, Eq.~(4.35)
	jacTmp[1][0] = -c*(*dummy)[0];										//!< the offdiagonal element of the voltage-voltage block, Eq.~(4.35)

	jacTmp[0][1] = (jacISIexpS - jacISIexpM)/(1 - gamma);				//!< the diagonal element of the voltage-current block, Eq.~(4.35)
	jacTmp[1][1] = -c*(*dummy)[1];										//!< the offdiagonal element of the voltage-current block, Eq.~(4.35)

	jacTmp[0][2] = 0;													//!< the diagonal element of the current-voltage block, Eq.~(4.35)
	jacTmp[1][2] = gamma*c*(*dummy)[0];									//!< the offdiagonal element of the current-voltage block, Eq.~(4.35)

	jacTmp[0][3] = jacISIexpS;											//!< the diagonal element of the current-current block, Eq.~(4.35)
	jacTmp[1][3] = gamma*c*(*dummy)[1];									//!< the offdiagonal element of the current-current block, Eq.~(4.35)

	return jacTmp;
}

//! The diagonal elements of the single spike Jacobian for the neurons that are not postsynaptic to the spiking neuron.
vector<vector<reell> > LE_clif::calc_JacElem_self(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));

	jacTmp[0][0] = jacISIexpM;											//!< the diagonal element of the voltage-voltage block, Eq.~(4.35)
	jacTmp[0][1] = (jacISIexpS - jacISIexpM)/(1 - gamma);				//!< the diagonal element of the voltage-current block, Eq.~(4.35)
	jacTmp[0][2] = 0;													//!< the diagonal element of the current-voltage block, Eq.~(4.35)
	jacTmp[0][3] = jacISIexpS;											//!< the diagonal element of the current-current block, Eq.~(4.35)

	return jacTmp;
}

//! The diagonal elements of the single spike Jacobian for the spiking neuron.
vector<vector<reell> > LE_clif::calc_JacElem_spiking(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));

	jacTmp[0][0] = jacISIexpM - (*dummy)[0]/tauM;									//!< the diagonal element of the voltage-voltage block, Eq.~(4.35)
	jacTmp[0][1] = (jacISIexpS - jacISIexpM)/(1 - gamma) - (*dummy)[1]/tauM;		//!< the diagonal element of the voltage-current block, Eq.~(4.35)
	jacTmp[0][2] = 0;																//!< the diagonal element of the current-voltage block, Eq.~(4.35)
	jacTmp[0][3] = jacISIexpS;														//!< the diagonal element of the current-current block, Eq.~(4.35)

	return jacTmp;
}

void LE_clif::dummy_Jacobian_spiking(vector<reell>* dummy)
{
	dummy->resize(2);

	reell tauMnu = tauM/(1 + Iext - state[0] + state[1]);				// 1-state[0] = 0

	(*dummy)[0] = -jacISIexpM*tauMnu;
	(*dummy)[1] = -(jacISIexpS - jacISIexpM)/(1 - gamma)*tauMnu;
}


//-------------------------- gamma = 2 --------------------------


LE_clif_fast2::LE_clif_fast2(reell tau) : LE_clif(tau)
{
	gamma = 2;
	tauS = tauM/2;
}

LE_clif_fast2::~LE_clif_fast2()
{

}

inline reell LE_clif_fast2::powGamma(reell x)
{
	return x*x;
}

inline void LE_clif_fast2::calc_spikeTime()
{
	if ((Iext <= 0) && (state[1] <= 0))								//!, this one is never gonna spike unless there is an excitatory input
		spikeTime = GSL_POSINF;
	else
	{
		//! this is the numerical method to calculate the roots of a square function from the numerical recipes book
		//! see thesis, p.~87
		reell a = -state[1];
		reell b = -(1 + Iext - state[0] - state[1]);
		reell c = Iext;

		reell d = -(b + GSL_SIGN(b)*sqrt(b*b - 4*a*c))/2;

		//! This defines two solutions d/a and c/d. Only one of them can be physically correct.
		reell x1 = c/d, x2 = d/a;

		//! Check if the logarithm of this root makes sense nd if the solution is unique.
		bool oneSol = false;

		if ((x1 > 0) && (x1 <= 1))
		{
			oneSol = true;
			spikeTime = -log(x1);
		}
		if ((x2 > 0) && (x2 <= 1))
		{
			if (oneSol)
			{
				cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
				throw(4);
			}
			else
			{
				oneSol = true;
				spikeTime = -log(x2);
			}
		}

		//! If both solutins are out of range, then the neuron will never spike until it gets another excitatory input, which
		//! is of course possible in the future evolution. For now, we'll set the spike time to infinity, though.
		if (!oneSol)
			spikeTime = GSL_POSINF;
	}

	calcSpikeTime = false;
}


//------------------------- gamma = 1/2 -------------------------


LE_clif_slow2::LE_clif_slow2(reell tau) : LE_clif(tau)
{
	gamma = 0.5;
	tauS = tauM*2;
}

LE_clif_slow2::~LE_clif_slow2()
{

}

inline reell LE_clif_slow2::powGamma(reell x)
{
	return sqrt(x);
}

inline void LE_clif_slow2::calc_spikeTime()
{
	if ((Iext <= 0) && (state[1] <= 0))								//!, this one is never gonna spike unless there is an excitatory input
		spikeTime = GSL_POSINF;
	else
	{
		//! this is the numerical method to calculate the roots of a square function from the numerical recipes book
		//! see thesis, p.~87
		reell a = -(1 + Iext - state[0] + 2*state[1]);
		reell b = 2*state[1];
		reell c = Iext;

		reell d = -(b + GSL_SIGN(b)*sqrt(b*b - 4*a*c))/2;

		//! This defines two solutions d/a and c/d. Only one of them can be physically correct.
		reell x1 = c/d, x2 = d/a;

		//! Check if the logarithm of this root makes sense nd if the solution is unique.
		bool oneSol = false;

		if ((x1 > 0) && (x1 <= 1))
		{
			oneSol = true;
			spikeTime = -2*log(x1);
		}
		if ((x2 > 0) && (x2 <= 1))
		{
			if (oneSol)
			{
				cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
				throw(4);
			}
			else
			{
				oneSol = true;
				spikeTime = -2*log(x2);
			}
		}

		//! If both solutins are out of range, then the neuron will never spike until it gets another excitatory input, which
		//! is of course possible in the future evolution. For now, we'll set the spike time to infinity, though.
		if (!oneSol)
			spikeTime = GSL_POSINF;
	}

	calcSpikeTime = false;
}


//-------------------------- gamma = 3 --------------------------


LE_clif_fast3::LE_clif_fast3(reell tau) : LE_clif(tau)
{
	gamma = 3;
	tauS = tauM/3;

	//! we need the machine precision for the calculation of the spike time
	MACHINE_PRECISION = numeric_limits<reell>::epsilon();
}

LE_clif_fast3::~LE_clif_fast3()
{

}

inline reell LE_clif_fast3::powGamma(reell x)
{
	return x*x*x;
}

inline void LE_clif_fast3::calc_spikeTime()
{
	if ((Iext <= 0) && (state[1] <= 0))								//!, this one is never gonna spike unless there is an excitatory input
		spikeTime = GSL_POSINF;
	else
	{
		//! this is the numerical method to calculate the roots of a cubic function from the numerical recipes book
		reell a = -state[1]/2;
		reell b = -(1 + Iext - state[0] - state[1]/2);
		reell c = Iext;

		if (abs(a) <= MACHINE_PRECISION)
		{
			reell x0 = -c/b;
			if ((x0 > 0) && (x0 <= 1))
				spikeTime = -log(x0);
			else
			{
				cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
				throw(4);
			}
		}
		else
		{
			reell Q = -b/a/3;
			reell R = c/a/2;

			reell R2 = R*R;
			reell Q3 = Q*Q*Q;

			if (R2 < Q3)
			{
				//there are 3 solutions, hopefully only one lies between 0 and 1 :)

				reell t = acos(R/sqrt(Q3));

				reell x1 = -2*sqrt(Q)*cos(t/3);
				reell x2 = -2*sqrt(Q)*cos((t + 2*M_PI)/3);
				reell x3 = -2*sqrt(Q)*cos((t - 2*M_PI)/3);

				bool oneSol = false;

				if ((x1 > 0) && (x1 <= 1))
				{
					oneSol = true;
					spikeTime = -log(x1);
				}
				if ((x2 > 0) && (x2 <= 1))
				{
					if (oneSol)
					{
						cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
						throw(4);
					}
					else
					{
						oneSol = true;
						spikeTime = -log(x2);
					}
				}
				if ((x3 > 0) && (x3 <= 1))
				{
					if (oneSol)
					{
						cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
						throw(4);
					}
					else
					{
						oneSol = true;
						spikeTime = -log(x3);
					}
				}
			}
			else
			{
				reell A = (R > 0) ? -cbrt(R + sqrt(R2-Q3)) : cbrt(-R + sqrt(R2-Q3));

				reell x4 = A + Q/A;

				if ((x4 > 0) && (x4 <= 1))
					spikeTime = -log(x4);
				else
					spikeTime = GSL_POSINF;
			}
		}
	}

	calcSpikeTime = false;
}


//------------------------- gamma = 1/3 -------------------------


LE_clif_slow3::LE_clif_slow3(reell tau) : LE_clif(tau)
{
	gamma = 1./3;
	tauS = tauM*3;

	//! we need the machine precision for the calculation of the spike time
	MACHINE_PRECISION = numeric_limits<reell>::epsilon();
}

LE_clif_slow3::~LE_clif_slow3()
{

}

inline reell LE_clif_slow3::powGamma(reell x)
{
	return cbrt(x);
}

inline void LE_clif_slow3::calc_spikeTime()
{
	if ((Iext <= 0) && (state[1] <= 0))								//!, this one is never gonna spike unless there is an excitatory input
		spikeTime = GSL_POSINF;
	else
	{
		//! this is the numerical method to calculate the roots of a cubic function from the numerical recipes book
		reell a = -(1 + Iext - state[0] + 3*state[1]/2);
		reell b = 3*state[1]/2;
		reell c = Iext;

		if (abs(a) <= MACHINE_PRECISION)
		{
			reell x0 = -c/b;
			if ((x0 > 0) && (x0 <= 1))
				spikeTime = -log(x0);
			else
			{
				cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
				throw(4);
			}
		}
		else
		{
			reell Q = -b/a/3;
			reell R = c/a/2;

			reell R2 = R*R;
			reell Q3 = Q*Q*Q;

			if (R2 < Q3)
			{
				//there are 3 solutions, hopefully only one lies between 0 and 1 :)

				reell t = acos(R/sqrt(Q3));

				reell x1 = -2*sqrt(Q)*cos(t/3);
				reell x2 = -2*sqrt(Q)*cos((t + 2*M_PI)/3);
				reell x3 = -2*sqrt(Q)*cos((t - 2*M_PI)/3);

				//cout << x1 << " " << x2 << " " << x3<< endl;

				bool oneSol = false;

				if ((x1 > 0) && (x1 <= 1))
				{
					oneSol = true;
					spikeTime = -3*log(x1);
				}
				if ((x2 > 0) && (x2 <= 1))
				{
					if (oneSol)
					{
						cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
						throw(4);
					}
					else
					{
						oneSol = true;
						spikeTime = -3*log(x2);
					}
				}
				if ((x3 > 0) && (x3 <= 1))
				{
					if (oneSol)
					{
						cout << "couldn't find a unique solution for this neuron with voltage=" << state[0] << " and current=" << state[1] << endl;
						throw(4);
					}
					else
					{
						oneSol = true;
						spikeTime = -3*log(x3);
					}
				}
			}
			else
			{
				reell A = (R > 0) ? -cbrt(R + sqrt(R2-Q3)) : cbrt(-R + sqrt(R2-Q3));

				reell x4 = A + Q/A;

				if ((x4 > 0) && (x4 <= 1))
					spikeTime = -3*log(x4);
				else
					spikeTime = GSL_POSINF;
			}
		}
	}

	calcSpikeTime = false;
}

