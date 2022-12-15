/*
 * LErapidtheta.h
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#ifndef LERAPIDTHETA_H_
#define LERAPIDTHETA_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>

using namespace std;

//! The rapid theta neuron class.

class LE_rapidtheta : public LE_phaseneuron
{

	public:
		LE_rapidtheta(reell, reell);
		LE_rapidtheta(reell, reell, reell, reell);
		virtual ~LE_rapidtheta();

		virtual void set_externalCurrent(reell&);

	protected:
		inline virtual reell PTC(reell);							//!< phase transition curve
		inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve

		reell rap;												//!< action potential onset rapidness
		reell AS;
		reell phiG;
};

#endif /* LERAPIDTHETA_H_ */
