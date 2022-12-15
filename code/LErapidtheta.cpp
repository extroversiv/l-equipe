/*
 * LErapidTheta.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#include "LErapidtheta.h"

LE_rapidtheta::LE_rapidtheta(reell r, reell tau) : LE_phaseneuron(tau, vector<reell> (1, -M_PI), vector<reell> (1, M_PI))
{
	//! This is the rapid theta neuron model. For r=1 it is the standard theta neuron model.

	//! Set the AP onset rapidness
	rap = r;

	//! The stable and the unstable fixed point are at -0.5 and +0.5 in the dimensionless voltage representation, respectively.
	Irheo = rap/(rap + 1)/2;					//!< rheo base current, Eq.~(3.3)
	AS = (rap + 1)/rap/2;						//!< curvature in stable area, Eq.~(3.4)
	phiG = M_PI*(rap - 1)/(rap + 1);			//!< glue point in phase description, Eq.~(3.13)

	//cout << "rap: " << rap << "\t Irheo: " << Irheo << "\tas: " << AS << "\tphiG: " << phiG << endl;
}


LE_rapidtheta::LE_rapidtheta(reell r, reell resetPhase, reell thresholdPhase, reell tau) : LE_phaseneuron(tau, vector<reell> (1, resetPhase), vector<reell> (1, thresholdPhase))
{
	//! This version of the rapid neuron model uses custom reset and threshold values. It can be thought of the rQIF, rapid
	//! quadratic integrate and fire. For r=1 its the standard QIF. For large resets it behaves like the rapid theta neuron model.

	//! The rest from here is a copy of the rapidtheta neuron constructor.
	//! Set the AP onset rapidness
	rap = r;

	//! The stable and the unstable fixed point are at -0.5 and +0.5 in the dimensionless voltage representation, respectively.
	Irheo = rap/(rap + 1)/2;					//!< rheo base current, Eq.~(3.3)
	AS = (rap + 1)/rap/2;						//!< curvature in stable area, Eq.~(3.4)
	phiG = M_PI*(rap - 1)/(rap + 1);			//!< glue point in phase description, Eq.~(3.13)

	//cout << "rap: " << rap << "\t Irheo: " << Irheo << "\tas: " << AS << "\tphiG: " << phiG << endl;
}


LE_rapidtheta::~LE_rapidtheta()
{

}

//! The external current sets the phase velocity.
void LE_rapidtheta::set_externalCurrent(reell& Iext)
{
	this->Iext = Iext;
	w = 2*sqrt(Iext/AS);						//!< Eq.~(3.12)
	calcSpikeTime = true;
}

//! The phase transition curve of the rapid theta neuron, Eq.~(3.16, 3.17).
inline reell LE_rapidtheta::PTC(reell c)
{
	reell phi = (state[0] - phiG)/2;
	c /= w/2;

	if (phi <= 0)
	{
		reell tanPhi = tan(AS*phi);

		if ((c > 0) && (-tanPhi < c))							//jumps from stable into unstable area
			phi = 2*atan(rap*(tanPhi + c))/AS/rap;
		else 													//stays in stable area
			phi = 2*atan(tanPhi + c)/AS;
	}
	else
	{
		reell tanPhi = tan(rap*AS*phi);

		if ((c <= 0) && (tanPhi < -rap*c))
			phi = 2*atan(tanPhi/rap + c)/AS;
		else
			phi = 2*atan(tanPhi + rap*c)/AS/rap;
	}

	return phi + phiG;	//new state
}

//! The derivative of the phase transition curve of the rapid theta neuron, Eq.~(3.19, 3.20).
inline reell LE_rapidtheta::PTCprime(reell c)
{
	reell phi = (state[0] - phiG)/2;
	c /= w/2;

	if (phi <= 0)
	{
		reell tanPhi = tan(AS*phi);

		reell den_part = 0;										//denominator part in derivative of phase transition curve

		if ((c > 0) && (-tanPhi < c))							//jumps from stable into unstable area
			den_part = rap*(tanPhi + c);
		else 													//stays in stable area
			den_part = tanPhi + c;

		return (tanPhi*tanPhi + 1)/(den_part*den_part + 1);
	}
	else
	{
		reell tanPhi = tan(rap*AS*phi);

		reell den_part = 0;

		if ((c <= 0) && (tanPhi < -rap*c))
			den_part = tanPhi/rap + c;
		else
			den_part = tanPhi + rap*c;

		return (tanPhi*tanPhi + 1)/(den_part*den_part + 1);
	}
}

