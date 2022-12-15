/*
 * LElif.h
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#ifndef LELIF_H_
#define LELIF_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>

using namespace std;

//! The leaky integrate and fire class

class LE_lif : public LE_phaseneuron
{
	public:
		LE_lif(reell);
		virtual ~LE_lif();

		virtual void set_externalCurrent(reell&);

	protected:
		inline virtual reell PTC(reell);							//!< phase transition curve
		inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve
};

#endif /* LELIF_H_ */
