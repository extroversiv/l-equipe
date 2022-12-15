/*
 * LEtwoDlinear.h
 *
 *  Created on: 14.02.2012
 *      Author: mik
 */

#ifndef LETWODLINEAR_H_
#define LETWODLINEAR_H_

#include "LEneuron.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "structures.h"
#include <limits>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;

//! The 2D linear neuron class.

class LE_twoDlinear : public LE_neuron
{
	public:
		LE_twoDlinear(reell, st_twoDlinear);
		virtual ~LE_twoDlinear();

		virtual reell get_phase();

		virtual void set_phase(reell&);
		virtual void set_externalCurrent(reell&);

		virtual void dummy_evolve(reell, vector<reell>*);
		virtual void evolve_dt(reell dt, vector< reell >*);
		virtual void evolve_spike(reell);
		virtual void reset();
		virtual void dummy_PostUpdate();

		virtual vector<vector<reell> > calc_JacElem_postsynaptic(vector<reell>*, reell);
		virtual vector<vector<reell> > calc_JacElem_self(vector<reell>*);
		virtual vector<vector<reell> > calc_JacElem_spiking(vector<reell>*);
		virtual void dummy_Jacobian_spiking(vector<reell>*);


		// gsl root finding
		static void proxyCompfdf(double tval, void *rootparams, double *y, double *dy);
		static void proxyRealfdf (double tval, void *rootparams, double *y, double *dy);
		static double proxyCompf  (double tval, void *rootparams);
		static double proxyRealf(double tval, void *rootparams);
		static double proxyCompdf (double tval, void *rootparams);
		static double proxyRealdf  (double tval, void *rootparams);

	protected:
	  
		//for home made root finding
		struct st_ftuple  //when using home made newton method
		{
			reell fval;
			reell dfval;
		};
		st_ftuple fdfval; 								//!>holds the root function value and its derivative
		typedef void (LE_twoDlinear::*HomemadeRtFn)(st_ftuple&,reell,reell,reell);  	//!pointer-to-function type for FReal/Fcomplex
		void FReal(st_ftuple&,reell,reell,reell);
		void FComplex(st_ftuple&,reell,reell,reell); 

		
		//for gsl root finding
		struct RootParams  //!>holds the root fuinction parameters
		{
		  double C1, C2; 
		  LE_twoDlinear* LE2DlinearOBJ;
		};
		RootParams rootparams;
		void rootRealfdf (double tval, void *rootparams, double *y, double *dy);
		void rootComplexfdf (double tval, void *rootparams, double *y, double *dy);
		double rootRealf (double tval, void *rootparams);
		double rootRealdf (double tval, void *rootparams);
		double rootComplexf (double tval, void *rootparams);
		double rootComplexdf (double tval, void *rootparams);

		

		
		virtual void calc_spikeTime();
		
		reell MACHINE_PRECISION;
		
		//!>renamed model paras
		
		reell alp;
		reell bet;
		reell deltauS;
		reell gam;
		reell Cw;		//!>steady state current ()reell Cv called Iext thoughout)			
		reell Iext;		//!< the external current																								
		
		//!internal model paras
		
		reell realPart; 		//!>1st term of eigen value
		reell omegaMOD; 	//!>second term of eigen value
		bool iscomplex; 	//!>flags case of complex eigen values 
		reell lambdaPlus; 	//!>positive conjugate eigen value
		reell lambdaMinus;	//!>negative conjugate eigen value


		reell kappaV;		//!>voltage fixed point
		reell kappaW;		//!>current fixed point
		reell A[4]; 		//!>differential matrix, used in propspike and storing jacobian biulding blocks, store as vector A=(A11,A12,A21,A22)							
 		reell dzdt[4];		//!>the difference in state velocities before and after spike

 		reell C3;		//!>used in constructor and calcspikeTime,just function of model paras
		reell Cv;		//!>steady state voltage (external current)

		reell initSpikingV;	//!>stores the initial state of the spikign neuron
		reell initSpikingW;	//!>stores the initial state of the spikign neuron
		reell jacISI;	//!>transmits spike time from eolvedtdummy to jacobiandummy
		
};

#endif /* LETWODLINEAR_H_ */
