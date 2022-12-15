/*
 * LEtype1type2.cpp
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#include "LEtype1type2.h"

LE_type1type2::LE_type1type2(reell a, reell tau) : LE_phaseneuron(tau, vector<reell> (1, -M_PI), vector<reell> (1, M_PI))
{
	/** This class interpoltes the type I and the type II PRCs by changing a
	 *  a = 0   -> PRC = 1 + cos(phi) -> type I
	 *  a = 0.5 -> PRC = sin(phi)     -> type II
	 *
	 *  The PRC is here $Z(\phi) = \cos(\phi - a\pi) + \cos(a\pi)$
	 */

	if ((a < 0) || (a > 0.5))
	{
		cout << "a = " << a;
		cout << ", but a needs to be between 0 and 0.5!" << endl;
		throw(1);
	}

	aPi = a*M_PI;
	aPiCos = cos(a*M_PI);

}

LE_type1type2::~LE_type1type2()
{

}

//! The external current sets the phase velocity.
void LE_type1type2::set_externalCurrent(reell& Iext)
{
	w = 2*sqrt(Iext);
	calcSpikeTime = true;
}

inline reell LE_type1type2::PTC(reell c)
{
	c /= w/2;

	return state[0] + c*(cos(state[0] - aPi) + aPiCos);
}

inline reell LE_type1type2::PTCprime(reell c)
{
	c /= w/2;

	return 1 - c*sin(state[0] - aPi);
}
