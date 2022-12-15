/*
 * LEtwoDlinear.cpp
 *
 *  Created on: 14.02.2012
 *      Author: max & mik
 */

#include "LEtwoDlinear.h"

LE_twoDlinear::LE_twoDlinear(reell tauM, st_twoDlinear paras) : LE_neuron(tauM, vector<reell> (2, 0), vector<reell> (2, 1))
{
	//! we need the machine precision for the calculation of the spike time
	MACHINE_PRECISION = numeric_limits<reell>::epsilon();

	//! The rheobase current is always 1 in the LIF model.
	//! The reset value is 0 and the threshold is 1.
	Irheo = 1;

	//! The state variable is two-dimensional in the twoDlinear model.
	state = vector<reell> (2);
	
	/**
	 * recall:	alpha	beta	gamma	delta	tauS	Cw
	 * cLIF		 1	0	0	1	0.5	0
	 * RF		-1	1	1	0	2	(0,1)
	 */



	if (paras.tauS == tauM)
	{
		cout << "tauS and tauM cannot be identical." << endl;
		throw(1);
	}

	reell tauS = paras.tauS/tauM;				//!< tauS in units of tauM

	alp    	=  paras.alpha;
	bet    	=  paras.beta;
	gam    	=  paras.gamma;
	deltauS	=  paras.delta/tauS;
	Cw	=  paras.Cw;
	



	//! global variable for the calculation of the real and complex eigenvalues

	//real part
	realPart = -(1 + 1/tauS)/2;

	//complex part
	reell radicandOmega	= realPart*realPart - (1 - alp*bet)/tauS;
	omegaMOD = sqrt(fabs(radicandOmega));
	iscomplex = (radicandOmega < 0) ? true : false;

	//eigenvalues
	lambdaPlus = realPart + omegaMOD;
	lambdaMinus = realPart - omegaMOD;
		



	//!DIfferential Matrix:used in propspike and storing jacobian biulding blocks, store as vector A=(A11,A12,A21,A22)
	A[0] = -1; 							//!>dissipativve voltage
	A[1] = alp; 						//!>current coupling to voltage
	A[2] = bet/tauS;					//!?voltage couplin gto current
	A[3] = -1/tauS;						//!>disipative current


	//!initialize the 3 biulding blocks of jacobian:
	//!1)the 2x2 linear operator, S(deltat), will be define das dummy
	//!2)two 2x1 vectors, dzdt, 2a) one for spiking and 2b) one for post syn pop, but the latter without J since J is neuron specific
	dzdt[0] = stateThreshold[0]*(A[0] + A[1]*bet/tauS);
	dzdt[1] = stateThreshold[0]*(A[2] + A[3]*bet/tauS);
	dzdt[2] = -(A[0]*gam + A[1]*deltauS);//!withoutJ, since its value depends on neuron...but will only use on post syn pop so doesn't matter...
	dzdt[3] = -(A[2]*gam + A[3]*deltauS);//!withoutJ,
	//!3) dtaudz, will be defined as dummy 
}

LE_twoDlinear::~LE_twoDlinear()
{

}

void LE_twoDlinear::set_phase(reell&)
{
	//state[0]=f(omega*phase)  f,g obtained by solving 2 equations that are linear in the states
	//state[1]=g(state[0])
	cout << "reell LE_twoDlinear::set_phase() not implemented yet!" << endl;
	throw(2);
	
}

reell LE_twoDlinear::get_phase()
{	    
	    //to=calc_ResetTime();
	    //t1=calc_SpikeTime();
  
	    //omega=1/(t1-to);
	    //phase=-omega*to;
  
	cout << "reell LE_twoDlinear::get_phase() not implemented yet!" << endl;
	throw(2);
}

void LE_twoDlinear::set_externalCurrent(reell& Iext)
{
	this->Iext = Iext;
	
	//!once Iext is defined we can calculate all remaining model parameters.
	
	Cv = Irheo + Iext;

	kappaV = (Cv + alp*Cw)/(1 - alp*bet);
	kappaW = (bet*Cv + Cw)/(1 - alp*bet);

	C3 = (kappaV - stateThreshold[0]);
	
	//if C3 less than 0, outside of restrictive conditions for root finding with tguess
	if (C3 < 0)
	{
	    cout << "C3 = " << C3 << " < 0 ! Root finding will fail" << endl;
	    throw(1);
	}

	if (iscomplex==false)
		C3 *= 2; //!just the way C3 is defined in the real case
	
	calcSpikeTime = true;
}

inline void LE_twoDlinear::calc_spikeTime()
{
	double tol = MACHINE_PRECISION;	//!>function value tolerance
	
	reell C1,C2; 
	if (state[0] < 1)
	{
		C1 = state[0] - kappaV;
		C2 = (alp*(state[1] - kappaW) - (1 + realPart)*C1)/omegaMOD;
		
		//!assign function parameters and define root function object
// 		rootparams.C1=C1;
// 		rootparams.C2=C2;
// 		rootparams.LE2DlinearOBJ=this;
// 		gsl_function_fdf FDF;
// 		FDF.params = (void*) &rootparams;

		
		//!The following assigns the correct (real/complex))root function, and 
		//!computes tguess to fall in interval in which points converge to correct root
		
		HomemadeRtFn fPoint;
		double tguess,tlower,tupper;
		bool bisectFlag=false;			//! flagwhether or not to run bisection, iinstead of Newton

		
		if (iscomplex == true)		//for complex eigen values
		{
			//commment out for homemade
// 			FDF.fdf = &proxyCompfdf;
// 			FDF.f = &proxyCompf;
// 			FDF.df = &proxyCompdf;

			fPoint = &LE_twoDlinear::FComplex;
			
			
			//!find extrema between -pi/2 and pi/2
			double temp1 = omegaMOD*C2 + realPart*C1;
			double temp2 = realPart*C2 - omegaMOD*C1;
			double tExtrema = atan( -temp1/temp2 )/omegaMOD;
			
			//!if tExtrema negative then look at next extrema time (which is necessarily positive)
			tExtrema += (tExtrema < 0) ? M_PI/omegaMOD : 0;
			
			
			if ((C2 > 0) || (realPart*C1 > -omegaMOD*C2))		//tExtrema is a maximum, 
			{
				tlower = 0;
				tupper = tExtrema;
				
				//inflection point
				double t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
				//ensure positive
				t_inflect += (t_inflect < 0)? M_PI/omegaMOD : 0;
				
				//is first positive inflection point in this region?
				if ( (tlower < t_inflect) && (t_inflect < tupper) )
				{
					(this->*fPoint)(fdfval, t_inflect, C1, C2);
					double f_inflect=fdfval.fval;
					
					if (f_inflect < -100*tol) 	//place to the right of inflection point
					{
						tlower = t_inflect;
						tguess = tlower+10*tol;
					}
					else if (f_inflect > 100*tol)	//place to the left of inflection point
					{
						tupper = t_inflect;
						tguess = tupper-10*tol;
					}
					else //root is near inflection point so derivative methods will fail, run bisection
					{
						bisectFlag = true;
						cout<<"BISECT!!"<<endl;
					}
				}
				else //no inflection point between 0 and tExtrema
				{
					tguess = tlower;
				}
			}
			else							//tExtream is a minimum
			{
				//valid region is between tExtrema=tmin and tExtrema+pi/omega
				tlower = tExtrema; //min
				tupper = tExtrema + M_PI/omegaMOD;//max
				
				//there must be an inflection point in this region
				double t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
				t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
				
				(this->*fPoint)(fdfval, t_inflect, C1, C2);
				double f_inflect=fdfval.fval;
				
				if (f_inflect < -100*tol)
				{
					tlower = t_inflect;
					tguess = tlower+10*tol;
				}
				else if (f_inflect > 100*tol)
				{
					tupper = t_inflect;
					tguess = tupper-10*tol;
				}
				else //root is near inflection point so derivative methods will fail, run bisection
				{
					bisectFlag = true;
					cout<<"BISECT!!"<<endl;
				}
			}
		}
		else //real eigen values
		{
			//comment out for homemade
// 			FDF.fdf=&proxyRealfdf;
// 			FDF.f=&proxyRealf;
// 			FDF.df=&proxyRealdf;

			fPoint = &LE_twoDlinear::FReal;	
			
			if (C2<C1) 		//!the extrema exists and is a min
			{
				//!place tguess between inflection point and root, where function is mono-convex/concave so derivative-based methods work
				double temp = lambdaMinus/lambdaPlus;
				double t_inflect	= log(temp*temp * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
				double tmin 		= log(   temp   * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
				
				(this->*fPoint)(fdfval, t_inflect, C1, C2);
				double f_inflect = fdfval.fval;
				
				if (f_inflect < -100*tol)
				{
					tlower = (t_inflect<0)? 0:t_inflect;
					tupper = GSL_POSINF;//potentially creates problems when swithcing to bisection
					tguess = tlower;
				}
				else if (f_inflect > 100*tol)
				{
					tlower = (tmin<0)? 0 : tmin;
					tupper = t_inflect;
					tguess = tupper-10*tol;
				}
				else //!root is near inflection point so derivative methods will fail, run bisection
				{
					bisectFlag = true;
					throw(1);
				}
// 				cout<<"flag1"<< endl;
			}
			else if (C2 > -C1) 	//! the extrema exists and is a max
			{
				
				//! max must occur at positive times since value at t=0 negative and asymptote is positive
				//! max must be positive, so region between t=0 and root is mono-concave

				tlower = 0;
				tupper = log( lambdaMinus/lambdaPlus*(C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
// 				cout<<"flag 2: "<< tupper << " "<<tlower <<endl;
				tguess = tlower;

			}
			else  			//!there is no extrema, then function monotonic so any value will do
			{
				tlower = 0;
				tupper = GSL_POSINF;
// 				cout<<"flag 3"<<endl;
				
				tguess = 10*tol;
			}
		}
		
// 		//!uncomment to initialize gsl algorithm
// 		const gsl_root_fdfsolver_type *T;
// 		gsl_root_fdfsolver *s;
// 		T = gsl_root_fdfsolver_newton;  //also try 
// // 		T = gsl_root_fdfsolver_steffenson;
// 		s = gsl_root_fdfsolver_alloc (T);
// 		gsl_root_fdfsolver_set (s, &FDF, tguess);
// 		printf ("using %s method\n", gsl_root_fdfsolver_name (s));
// 		int statusRootval;
// 		int statusFnval;

		//!run  algorithm
		int iter = 0, max_iter = 100;
		double told, tnew = tguess;//,fval;
		double toldold;
		//!run through algorithm
		cout.precision(18);
		
		do
		{
			if (tupper<tlower)
			{
				cout<< iter<<": brackets inverted!"<<endl;
				cout<< C1 << " "<<C2<<" "<<C3<<endl;
				throw(1);
			}
			if ( (tnew < tlower) || (tnew > tupper) )
			{
				cout<< iter<<": tnew out of bracket"<<endl;
				throw(1);
			}	
			
			  
			iter++;
			
			//!if root doesn't converge by 10 iterations (usually because oscilalting aroudn precision), switch to bisection
			if (iter > 10) //|| (fabs(told-tnew)>10*tol)
				bisectFlag = true;	

			toldold = told;
			//compute function value used to get next root estimate
			told = tnew;
			(this->*fPoint)(fdfval, told, C1, C2);
			
			//update brackets
			tlower = (fdfval.fval>=0)? tlower:told;
			tupper = (fdfval.fval>=0)? told:tupper;

			//homemade newton:
			if (bisectFlag == false)
			{
				tnew = told - (fdfval.fval)/(fdfval.dfval);
				
				//gsl:
	// 			fval=rootRealf(tnew, &rootparams);
	// 			statusRootval = gsl_root_fdfsolver_iterate (s);
	// 			tnew = gsl_root_fdfsolver_root (s);
	// 			statusRootval = gsl_root_test_delta (tnew, told, tol, 0); //epsabs, epsrel: |x_1 - x_0| < epsabs + epsrel |x_1|
	// 			statusFnval = gsl_root_test_residual(fval,10*tol);	
				
				//! apply bisection if tnew goes out of bounds.
				if ( (tnew < tlower) || (tnew > tupper) )
				{
					tnew = (tupper + tlower)/2;
					
// 					cout<<"bisect in Newton"<<endl;
// 					printf ("%5d %10.18f %+10.18f %10.18f\n", iter, tupper-tlower, tnew-told,fdfval.fval);
				}
			}
			else
			{
				//correct for case that algorithm only approached root in +/- direction, so upper/lower bound never reassigned
				//root approached from below
				if (tupper == GSL_POSINF)
				{
// 					cout<<"correct for posinf "<<endl;
					tupper = toldold;
					do 
					{
						//guess a close upper bound
						tupper += 10*fabs(told-toldold);
						(this->*fPoint)(fdfval, tupper, C1, C2);
					}			
					while ( fdfval.fval < 0 ); 
				}
					
				
				//! apply bisection
 				
				tnew= (tupper + tlower)/2;
				
// 				cout<<"bisect at end"<<endl;
// 				printf ("%5d %10.18f %+10.18f %10.18f\n", iter, tupper-tlower, tnew-told,fdfval.fval);

			}	
		}	
		while ( ( (fabs(tnew-told)> tol) && (fabs(tlower-tupper)>tol) ) && (iter <= max_iter) );		//! homemade conditions && (fabs(fdfval.fval)> tol) 
		
		spikeTime=tnew;												//! homemade estiamte

// 		while ( (statusFnval == GSL_CONTINUE) && (fabs(tupper-tlower) > tol) && (iter < max_iter);	//! gsl conditions
// 
// 		spikeTime= gsl_root_fdfsolver_root(s);								//! gsl estimate
// 		gsl_root_fdfsolver_free (s);									//! free gsl object
// 					
		if (!(iter < max_iter))
		{
			cout.precision(16);
			cout << " Maxed out iterations in root finding algorithm" << endl;
			cout << "state[0]" << state[0] << endl;
			cout << "state[1]" << state[1] << endl;
			cout << "C1 = " << C1 << endl;
			cout << "C2 = " << C2 << endl;
			cout << "C3 = " << C3 << endl;
			cout << "maschine precision="<<MACHINE_PRECISION << endl;
			throw(1);
		}
	}
	else if (state[0] > 1)
	{
		cout << "phase greater than 1! V = " << state[0] << " I = " << state[1] << endl;
		cout<<"Iext"<< this->Iext;
		throw(1);
	}
	else
	{
		spikeTime = 0;
	}
	
	if (spikeTime < -MACHINE_PRECISION)
	{
	    cout << "The neurons next spike time from root finding " << spikeTime << " was in the past!" << endl;
	    cout << "state[0]" << state[0] << endl;
	    cout << "state[1]" << state[1] << endl;
	    cout.precision(18);
	    cout << "C1=" << C1 << endl;
	    cout << "C2=" << C2 << endl;
	    cout << "C3=" << C3 << endl;
	    throw(1);
	}
	calcSpikeTime = false;
}

void LE_twoDlinear::dummy_evolve(reell dt, vector<reell>* dummy)
{
	//!the following calculates the propogation matrix,S, 
	//!it is only run for one (the 1st) neuron, but S is dimensionless and valid for all neurons (even with heterogeneous taus)
	//!The dummy variable is used in evolve_dt for all neurons.
	
	dt /= tauM;
	
	
	(*dummy).resize(4);

	reell S[4];

	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		S[0] = expRdt*( cosOmega + sinOmega*(A[0]-realPart)/omegaMOD );
		S[1] = expRdt*sinOmega*A[1]/omegaMOD;
		S[2] = expRdt*sinOmega*A[2]/omegaMOD;
		S[3] = expRdt*( cosOmega + sinOmega*(A[3]-realPart)/omegaMOD );
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		S[0] = coeff*A[0] + tempS;
		S[1] = coeff*A[1];
		S[2] = coeff*A[2];
		S[3] = coeff*A[3] + tempS;
	}

	//!assign to dummy
	(*dummy)[0] = S[0];
	(*dummy)[1] = S[1];
	(*dummy)[2] = S[2];
	(*dummy)[3] = S[3];

}

void LE_twoDlinear::evolve_dt(reell dt, vector<reell>* dummy)
{
	
	//! increments time between spikes in the network
	jacISI += dt;					
	
	//!evolve neurons by dt. dummy is the propogration matrix for this time.
	reell v_temp=state[0];
	state[0] = (*dummy)[0]*(state[0] - kappaV) + (*dummy)[1]*(state[1]-kappaW)+kappaV;
	state[1]  = (*dummy)[2]*(v_temp - kappaV)  + (*dummy)[3]*(state[1]-kappaW)+kappaW;

	//!adjust the neuron's spike time
	calcSpikeTime = true;
	
	//! We cannot just use spikeTime -= dt because this acccumulates precision losses. For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.
	

}

void LE_twoDlinear::evolve_spike(reell c)
{
	//! \a state[0] - the voltage and the current\a state[1] can change
	state[0] +=     gam*c;
	state[1] += deltauS*c;

	calcSpikeTime = true;
}

//! Set the state variable to the reset value.
void LE_twoDlinear::reset()
{
	
  
  
	//!voltage and current reset
	state[0]  = stateReset[0];
	state[1] -= stateReset[1];//reset value not clear, for resonator with slow currents, no rest is fine

	calcSpikeTime = true;
}

void LE_twoDlinear::dummy_PostUpdate()
{
	//! This is saved for each neuron, but only to be used in dummy_Jacobian_spiking
	initSpikingV = state[0];
	initSpikingW = state[1];
	
	//reset at network spike time
	jacISI=0;
  
}

void LE_twoDlinear::dummy_Jacobian_spiking(vector<reell>* dummy)
{
	//! the following calculates the propogation matrix and dtau/dz. It is run only for the spikign neuron.
	(*dummy).resize(6);

	reell dt = jacISI/tauM; //units of spiking neuron's 'membrane time constant
	
	reell dS11dt;
	reell dS12dt;
	reell S[4];
	
	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		S[0] = expRdt*( cosOmega + sinOmega*(A[0] - realPart)/omegaMOD );
		S[1] = expRdt*sinOmega*A[1]/omegaMOD;
		S[2] = expRdt*sinOmega*A[2]/omegaMOD;
		S[3] = expRdt*( cosOmega + sinOmega*(A[3] - realPart)/omegaMOD );
		
		dS11dt = realPart*S[0] - expRdt*(omegaMOD*sinOmega - cosOmega*(A[0] - realPart));
		dS12dt = realPart*S[1] + expRdt*cosOmega*A[1];
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		S[0] = coeff*A[0] + tempS;
		S[1] = coeff*A[1];
		S[2] = coeff*A[2];
		S[3] = coeff*A[3] + tempS;
		
		reell lambexpplus = lambdaPlus*expplus;
		reell lambexpminus = lambdaMinus*expminus;

		dS11dt = lambexpplus + (lambexpplus - lambexpminus)/(2*omegaMOD)*(A[0] - lambdaPlus);
		dS12dt = (lambexpplus - lambexpminus)*A[1]/(2*omegaMOD);
	}
	

// 	dS11dt /= tauM; //convert back to milliseconds as if it was the spiking neuron. assumes gamma same for all.
// 	dS12dt /= tauM;


	(*dummy)[0] = S[0];
	(*dummy)[1] = S[1];
	(*dummy)[2] = S[2];
	(*dummy)[3] = S[3];
	
	//dtaudz in temporal units of milliseconds
	(*dummy)[4] = S[0]/( -dS11dt/S[0]*(stateThreshold[0] - kappaV - S[1]*(initSpikingW - kappaW)) - dS12dt*(initSpikingW - kappaW) );
	(*dummy)[5] = S[1]/( -dS12dt/S[1]*(stateThreshold[0] - kappaV - S[0]*(initSpikingV - kappaV)) - dS11dt*(initSpikingV - kappaV) );
}

//! The diagonal elements of the single spike Jacobian for the postsynaptic neurons.
vector<vector<reell> > LE_twoDlinear::calc_JacElem_postsynaptic(vector<reell>* dummy, reell c)
{

	vector<vector<reell> > jacTmp(2, vector<reell>(4));				//!< the vector with Jacobian elements														
	//! 4 elements for the diagonal and the nondiagonal part, each.
	jacTmp[0][0] = (*dummy)[0];	
	jacTmp[0][1] = (*dummy)[1];
	jacTmp[0][2] = (*dummy)[2];		
	jacTmp[0][3] = (*dummy)[3];
	
// 	(*dummy)[4]/=tauM; //convert back to units of membrane time constant. This time of the post synaptic neuron.
// 	(*dummy)[5]/=tauM;
	
	jacTmp[1][0] = c*dzdt[2]*(*dummy)[4];															
	jacTmp[1][1] = c*dzdt[2]*(*dummy)[5];																						
	jacTmp[1][2] = c*dzdt[3]*(*dummy)[4];																		
	jacTmp[1][3] = c*dzdt[3]*(*dummy)[5];								
	return jacTmp;
}

vector<vector<reell> > LE_twoDlinear::calc_JacElem_self(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));
	//!the diagonal elements are just the propagation matrix, S
	jacTmp[0][0] = (*dummy)[0];											
	jacTmp[0][1] = (*dummy)[1];											
	jacTmp[0][2] = (*dummy)[2];											
	jacTmp[0][3] = (*dummy)[3];
	return jacTmp;
}

vector<vector<reell> > LE_twoDlinear::calc_JacElem_spiking(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));
	
// 	(*dummy)[4]/=tauM; //convert back to units of membrane time constant. This time of the spiking neuron.
// 	(*dummy)[5]/=tauM; 
	
	jacTmp[0][0] = (*dummy)[0] + dzdt[0]*(*dummy)[4];
	jacTmp[0][1] = (*dummy)[1] + dzdt[0]*(*dummy)[5];
	jacTmp[0][2] = (*dummy)[2] + dzdt[1]*(*dummy)[4];
	jacTmp[0][3] = (*dummy)[3] + dzdt[1]*(*dummy)[5];
	return jacTmp;
}



void LE_twoDlinear::FReal(LE_twoDlinear::st_ftuple &fdfval, reell tval, reell C1, reell C2)
{

//  	//! Some temporary variables
// 	reell fplus =  (C1 + C2)*exp(lambdaPlus*tval);
// 	reell fminus = (C1 - C2)*exp(lambdaMinus*tval);
// 
//  	//! The function value and the derivative.
// 	fdfval.fval  = fplus + fminus + C3;
// 	fdfval.dfval = fplus*lambdaPlus + fminus*lambdaMinus;
	
	reell expPlus  = exp(lambdaPlus*tval);
	reell expMinus = exp(lambdaMinus*tval);

	fdfval.fval  = C1*(expPlus + expMinus) + C2*(expPlus - expMinus) + C3;
	fdfval.dfval  = C1*(lambdaPlus*expPlus + lambdaMinus*expMinus) + C2*(lambdaPlus*expPlus - lambdaMinus*expMinus);

}

void LE_twoDlinear::FComplex(LE_twoDlinear::st_ftuple &fdfval, reell tval, reell C1, reell C2)
{
	//! Some temporary variables
	reell argO 		    = omegaMOD*tval;
	reell cosOmegatemp 	= cos(argO);
	reell sinOmegatemp 	= sin(argO);
	reell expRdttemp 	= exp(realPart*tval);
	reell temp 		    = expRdttemp*(C1*cosOmegatemp + C2*sinOmegatemp);
	
	//! The function value and the derivative.
	fdfval.fval = temp + C3;
	fdfval.dfval = realPart*temp - omegaMOD*expRdttemp*(C1*sinOmegatemp - C2*cosOmegatemp);
} 


void LE_twoDlinear::rootRealfdf (double tval, void *params, double *y, double *dy)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  double expPlus  = exp(lambdaPlus*tval);
  double expMinus = exp(lambdaMinus*tval);

  *y =  C1*(expPlus + expMinus) + C2*(expPlus - expMinus) + C3;
  *dy = C1*(lambdaPlus*expPlus + lambdaMinus*expMinus) + C2*(lambdaPlus*expPlus - lambdaMinus*expMinus);
}

void LE_twoDlinear::proxyRealfdf(double tval, void *params, double *y, double *dy)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	rootparams->LE2DlinearOBJ->rootRealfdf(tval, params, y, dy);
}


void LE_twoDlinear::rootComplexfdf (double tval, void *params, double *y, double *dy)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  //! Some temporary variables
  double argO 		    = omegaMOD*tval;
  double cosOmegatemp 	= cos(argO);
  double sinOmegatemp 	= sin(argO);
  double expRdttemp 	= exp(realPart*tval);
  double temp 		    = expRdttemp*(C1*cosOmegatemp + C2*sinOmegatemp);

  *y =  temp + C3;
  *dy = realPart*temp - omegaMOD*expRdttemp*(C1*sinOmegatemp - C2*cosOmegatemp);
}
void LE_twoDlinear::proxyCompfdf(double tval, void *params, double *y, double *dy)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	rootparams->LE2DlinearOBJ->rootComplexfdf(tval, params, y, dy);	
}


double LE_twoDlinear::rootRealf (double tval, void *params)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  double expPlus  = exp(lambdaPlus*tval);
  double expMinus = exp(lambdaMinus*tval);

  return C1*(expPlus + expMinus) + C2*(expPlus - expMinus) + C3;
  //*dy = C1*(lambdaPlus*expPlus + lambdaMinus*expMinus) + C2*(lambdaPlus*expPlus - lambdaMinus*expMinus);
}
double LE_twoDlinear::proxyRealf(double tval, void *params)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	return rootparams->LE2DlinearOBJ->rootRealf(tval, params);	
}


double LE_twoDlinear::rootComplexf (double tval, void *params)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  //! Some temporary variables
  double argO 		    = omegaMOD*tval;
  double cosOmegatemp 	= cos(argO);
  double sinOmegatemp 	= sin(argO);
  double expRdttemp 	= exp(realPart*tval);
  double temp 		    = expRdttemp*(C1*cosOmegatemp + C2*sinOmegatemp);

  return  temp + C3;
  //*dy = realPart*temp - omegaMOD*expRdttemp*(C1*sinOmegatemp - C2*cosOmegatemp);
}
double LE_twoDlinear::proxyCompf(double tval, void *params)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	return rootparams->LE2DlinearOBJ->rootComplexf(tval, params);	
}


double LE_twoDlinear::rootRealdf (double tval, void *params)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  double expPlus  = exp(lambdaPlus*tval);
  double expMinus = exp(lambdaMinus*tval);

  //*y =  C1*(expPlus + expMinus) + C2*(expPlus - expMinus) + C3;
  return C1*(lambdaPlus*expPlus + lambdaMinus*expMinus) + C2*(lambdaPlus*expPlus - lambdaMinus*expMinus);
}
double LE_twoDlinear::proxyRealdf(double tval, void *params)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	return rootparams->LE2DlinearOBJ->rootRealdf(tval, params);	
}


double LE_twoDlinear::rootComplexdf (double tval, void *params)
{
  struct RootParams *p = (struct RootParams *) params;

  double C1 = p->C1;
  double C2 = p->C2;

  //! Some temporary variables
  double argO 		    = omegaMOD*tval;
  double cosOmegatemp 	= cos(argO);
  double sinOmegatemp 	= sin(argO);
  double expRdttemp 	= exp(realPart*tval);
  double temp 		    = expRdttemp*(C1*cosOmegatemp + C2*sinOmegatemp);

  //*y =  temp + C3;
  return realPart*temp - omegaMOD*expRdttemp*(C1*sinOmegatemp - C2*cosOmegatemp);
}
double LE_twoDlinear::proxyCompdf(double tval, void *params)
{ 
	struct LE_twoDlinear::RootParams *rootparams = (struct LE_twoDlinear::RootParams *) params;
	return rootparams->LE2DlinearOBJ->rootComplexdf(tval, params);	
}





		//! gsl loop
// 		do 
// 		{
// 			iter++;
// 			statusRel = gsl_root_fdfsolver_iterate (s);
// 			tval0 = tval;
// 			tval = gsl_root_fdfsolver_root (s);
// 			
// 			
// 			if (tval < tmin)
// 			{tval = (tval0 + tmin)/2; cout <<"lowbound "<<tval<<endl;}
// 			if (tval > tmax)
// 			{tval = (tval0 + tmin)/2; cout <<"lowbound "<<tval<<endl;}
// 			
// 			statusRel = gsl_root_test_delta (tval, tval0, tol, 0); //epsabs, epsrel: |x_1 - x_0| < epsabs + epsrel |x_1|
// 			fval=rootRealf(tval, &rootparams);
// 			statusFn = gsl_root_test_residual(fval,10*tol);
// 			
// // 	      	  if ((statusRel == GSL_SUCCESS)& (statusFn == GSL_SUCCESS))
// // 	      	    printf ("Converged:\n");
// // 	      	  //fval=rootRealf(tval, &rootparams);
// // 	      	  printf ("%5d %10.18f %+10.19f %10.18f\n",
// // 	      		  iter, tval, tval - tval0,fval );
// 		      
// 		} 
// 		while (( !(statusFn == GSL_SUCCESS) && !(statusRel == GSL_SUCCESS) ) && (iter < max_iter)) ;
// 		spikeTime= gsl_root_fdfsolver_root(s);
// 		gsl_root_fdfsolver_free (s);//does this run if it is after the return?