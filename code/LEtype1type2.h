/*
 * LEtype1type2.h
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#ifndef LETYPE1TYPE2_H_
#define LETYPE1TYPE2_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>

using namespace std;

//! The phaseneuron class to interpolate freely between type I and type II phase response curve, see [Schleimer and Stemmler, PRL 103, 248105 (2009)].

class LE_type1type2 : public LE_phaseneuron
{
	public:
		LE_type1type2(reell, reell);
		virtual ~LE_type1type2();

		virtual void set_externalCurrent(reell&);

	protected:
		inline virtual reell PTC(reell);							//!< phase transition curve
		inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve

		reell aPi;
		reell aPiCos;
};

#endif /* LETYPE1TYPE2_H_ */
