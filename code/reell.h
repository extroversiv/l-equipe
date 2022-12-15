#ifndef REELL_H_
#define REELL_H_

#ifdef PAR
	#include "mpi.h"
#endif


#ifndef PAR

	//! redefine the type double (reelle Zahl) for an easy change to higher precision (typedef long double reell;)
	typedef double reell;

#else

	//! The typedefinition change is not working in the parallel version.
	//! Thus, don't change double to long double in the parallel version.
	typedef double reell; //Don't change this, got it?

#endif


#endif /*REELL_H_*/
