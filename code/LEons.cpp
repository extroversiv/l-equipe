/*
 * LEons.cpp
 *
 *  Created on: 09.01.2012
 *      Author: mik
 */

#include "LEons.h"

LE_ons::LE_ons(int col, int row, int dim, int seed) : COL(col), ROW(row), DIM(dim)
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


	//! initialize random number generator
	gsl_rng_env_setup();
	const gsl_rng_type *TYPE = gsl_rng_mt19937;
	rng = gsl_rng_alloc(TYPE);


	/** Standard orthonormal system (ONS) initiliazation ons[column][row*state dimension]
	 *  	column is equal to the number of vectors/Lyapunov exponents
	 *  	row is equal to the number of neurons on the loacl node (N/nP)
	 *  	state dimension is the dimension of the neuron's state variable
	 */


	// if (myID == 0) cout  << "cols: " << COL << "\trows: " << ROW << "\tdim: " << DIM << endl;

#ifdef BLAS
	// if (myID == 0) cout << "compiled with BLAS support for the orthogonalizations ... " << endl;
#endif

	//! ons initialized this way, to work correctly with BLAS!!!!!!
	//! This assures that the matrix is stored in one contiguous array
	ons = new reell* [COL];
	ons[0] = new reell [COL*ROW*DIM];
	for (int n=1; n<COL; n++)
		ons[n] = ons[n-1] + ROW*DIM;

	//! initialize the norm vector
	normONS = vector<reell> (COL);

	//! initialize the projections vectors used during the orthonormalization
	projections = vector<reell> (COL);

#ifdef PAR
	projectionsLocal = vector<reell> (COL);
#endif



	//! By default the projections should not be saved (This needs a lot of memory).
	saveProjections = false;


	//! Prepare projectionsStorage to hold the projections. This is a triangular matrix.
	if (myID == 0)
	{
		projStore = vector<vector<reell> > (COL);
		for (int i=0; i<COL; i++)
			projStore[i] = vector<reell> (COL-i);
	}


	//! Initialize the ons
	initONS(seed);

}

LE_ons::~LE_ons()
{
	gsl_rng_free(rng);

	delete[] ons[0];
	delete[] ons;

}

void LE_ons::initONS(int seed)
{
	gsl_rng_set(rng, seed);

	for (int n=0; n<COL; n++)
		for (int m=0; m<ROW*DIM; m++)
			ons[n][m] = gsl_rng_uniform(rng) - 0.5;				//!< positive and negative values, othonormalization afterwards

	orthonormalizeONS(COL);
}

void LE_ons::orthonormalizeONS(int Ncol)
{
	if (Ncol > COL)
	{
		cout << "The number of requested columns " << Ncol << " exceeds the number of orthonormal vectors " << COL << endl;
		throw(1);
	}

	/** There are many ways the QR othonormalization can be done.
	 *  For large matrices with parallel implementation, the row version of the modified Gram-Schmidt orhonormalization
	 *  turned out to be best. The parallel communication is low and efficient BLAS2 multiplications of the vectors can be used.
	 */

	//! effective number of rows
	int rowEff = ROW*DIM;


	/** BLAS speeds up the calculation by optimizing the array access to the cache of the processor.
	 *  Therefore, the specifically tuned package ATLAS should be used by linking with -lcblas -latlas.
	 *  An example is a calculation with N=500, K=50 of theta neurons for 7000 spikes and reorthonormalizations
	 *  took 1500 seconds without BLAS, 1300 s with BLAS provided by GSL and 1000 s with BLAS by ATLAS.
	 */

#ifdef BLAS

	//*************************
	//*      WITH BLAS        *
	//*************************
	// use -DBLAS compilerflag and link the library -lgslcblas for nonoptimzed BLAS or -lcblas -latlas for optimized BLAS support

	for(int i=0; i<Ncol; i++)
	{
		//! 1) calculate norm of vector i
		reell normI = cblas_dnrm2(rowEff, ons[i], 1);			//!< returns \sqrt(ons_trans*ons)

#ifdef PAR
		//! Add the parts of the norm from all other nodes.
		reell local, global;
		local = normI*normI;

		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		normI = sqrt(global);
#endif

		normONS[i] = normI;

		//! normalize vector i
		cblas_dscal(rowEff, 1/normI, ons[i], 1);				//!< scales ons[i] by 1/norm

		//! 2) othogonalize other vectors with i

		if (i < Ncol-1)
		{
			//!calculate projections
	#ifdef PAR
			cblas_dgemv(CblasRowMajor, CblasNoTrans, Ncol-(i+1), rowEff, 1, ons[i+1], rowEff, ons[i], 1, 0, &projectionsLocal[1], 1);	//!< returns projections = 1*ons[j>i]*ons[i] + 0*projections

			//! Add the parts from all other nodes.
			MPI_Allreduce(&projectionsLocal[1], &projections[1], Ncol-(i+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
			cblas_dgemv(CblasRowMajor, CblasNoTrans, Ncol-(i+1), rowEff, 1, ons[i+1], rowEff, ons[i], 1, 0, &projections[1], 1);		//!< returns projections = 1*ons[j>i]*ons[i] + 0*projections
	#endif

			//!subtract projections
			cblas_dger(CblasRowMajor, Ncol-(i+1), rowEff, -1, &projections[1], 1, ons[i], 1, ons[i+1], rowEff);							//!< returns ons[j>i] = -1*projections*ons[i] + ons[j>i]
		}

		//! Collect the projections of all vectors to store them in allSavedProjections for the covariant Lyapunov vector calculations
		if ((myID == 0) && saveProjections)
		{
			projections[0] = normI;							// this is the norm, sqrt(ons[i]*ons[i])

			//!store all projections(including the norm), in triangular form to save some memory
			projStore[i].assign(projections.begin(), projections.begin() + Ncol-i);
			// it's not Ncol-(i+1) here, because void assign ( InputIterator first, InputIterator last ); with the range [first,last)
			// so the last element is excluded in the assign method and the number of elements is thus Ncol-(i+1)
		}

	}

#else

	//*************************
	//*        NO BLAS        *
	//*************************

	for(int i=0; i<Ncol; i++)
	{
		//! calculate norm of vector i
		reell normI = 0;
		for(int k=0; k<rowEff; k++)
			normI += ons[i][k]*ons[i][k];

	#ifdef PAR
			//! Add the parts of the norm from all other nodes.
			reell global;
			MPI_Allreduce(&normI, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			normI = global;
	#endif
		normI = sqrt(normI);
		normONS[i] = normI;

		//! normalize vector i
		for(int k=0; k<rowEff; k++)
			ons[i][k] /= normI;

		//! othogonalize other vectors with i
	#ifdef PAR
		for(int j=i+1; j<Ncol; j++)
		{
			//! calculate projections
			projectionsLocal[j-i] = 0;
			for(int k=0; k<rowEff; k++)
				projectionsLocal[j-i] += ons[j][k]*ons[i][k];
		}

		//! Add the parts from all other nodes.
		MPI_Allreduce(&projectionsLocal[1], &projections[1], Ncol-(i+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
		for(int j=i+1; j<Ncol; j++)
		{
			//! calculate projections
			projections[j-i] = 0;
			for(int k=0; k<rowEff; k++)
				projections[j-i] += ons[j][k]*ons[i][k];
		}
	#endif

		//! subtract projections
		for(int j=i+1; j<Ncol; j++)
			for(int k=0; k<rowEff; k++)
				ons[j][k] -= projections[j-i]*ons[i][k];


		//! Collect the projections of all vectors to store them in allSavedProjections for the covariant Lyapunov vector calculations
		if ((myID == 0) && saveProjections)
		{
			projections[0] = normI;							// this is the norm, sqrt(ons[i]*ons[i])

			//!store all projections(including the norm), in triangular form to save some memory
			projStore[i].assign(projections.begin(), projections.begin() + Ncol-i);
			// it's not Ncol-(i+1) here, because void assign ( InputIterator first, InputIterator last ); with the range [first,last)
			// so the last element is excluded in the assign method and the number of elements is thus Ncol-(i+1)
		}


	}

#endif


	if ((myID == 0) && saveProjections)
		allSavedProjections.push_back(projStore);

}

void LE_ons::initCLV()
{
	if (myID == 0)
	{
		//! The covariant Lyapunov vectors are unique and must be independent of the initial seed, therefore seed is set to 1.
		gsl_rng_set(rng, 1);

		//! initialize the norm vector
		normCLV = vector<reell> (COL);

		//! initialize the coefficient of the covariant Lyapunov vectors w.r.t to the ONS
		clv = vector<vector<reell> > (COL);
		for (int i=0; i<COL; i++)
		{
			clv[i] = vector<reell> (i+1);

			for (int j=0; j<=i; j++)
				clv[i][j] = gsl_rng_uniform(rng) - 0.5;				//!< positive and negative values, othonormalization afterwards
		}

		normalizeCLV(COL);
	}
	else
	{
		cout << "The covariant Lyapunov vector calculation is solely done on root, but here it was called on processor #" << myID << endl;
		throw(1);
	}
}

void LE_ons::normalizeCLV(int Ncol)
{
	if (myID == 0)
	{

		if (Ncol > COL)
		{
			cout << "The number of requested columns " << Ncol << " exceeds the number of Lyapunov vectors " << COL << endl;
			throw(1);
		}

		for(int i=0; i<Ncol; i++)
		{
#ifdef BLAS
			//! calculate norm of vector i
			reell normI = cblas_dnrm2(i+1, &clv[i].front(), 1);				//!< returns \sqrt(clv.*clv)

			normCLV[i] = normI;

			//! normalize vector i
			cblas_dscal(i+1, 1/normI, &clv[i].front(), 1);
#else
			//! calculate norm of vector i
			reell normI = 0;
			for(int j=0; j<=i; j++)
				normI += clv[i][j]*clv[i][j];

			normI = sqrt(normI);
			normCLV[i] = normI;

			//! normalize vector i
			for(int j=0; j<=i; j++)
				clv[i][j] /= normI;
#endif
		}

	}
	else
	{
		cout << "The covariant Lyapunov vector calculation is solely done on root (sequential), but here it was called on processor #" << myID << endl;
		throw(1);
	}
}

void LE_ons::backwardIterationCLV(unsigned s, int Ncol)
{
	//! Calculates the backward iteration s, based on the saved projection matrix \aallSavedProjections.

	if (s >= allSavedProjections.size())
	{
		cout << "tried to calculate the backward iteration of step " << s << ", but there are only " << allSavedProjections.size() << " projections saved." << endl;
		throw(1);
	}

	if (s < 0)
	{
		cout << "tried to calculate the backward iteration of step " << s << ", this is wrong." << endl;
		throw(1);
	}

	//! 1) invert the projection matrix
	//! Gauss elimination by hand: start from the bottom right element and calculate inverse projection matrix (use \aprojectionsStorage to store the inverse matrix)
	for (int i=Ncol-1; i>=0; i--)
	{
		//! diagonal elements
		projStore[i][0] = 1/allSavedProjections[s][i][0];

		//! off-diagonal elements
		for (int j=i+1; j<Ncol; j++)
		{
			reell sum = 0;
			for (int k=i+1; k<=j; k++)
				sum += allSavedProjections[s][i][k-i] * projStore[k][j-k];		//no simple BLAS version because -k in last index

			projStore[i][j-i] = -sum/allSavedProjections[s][i][0];
		}
	}



	//! 2) multiply inverse projection matrix with coefficient matrix
	for (int j=0; j<Ncol; j++)
		for (int i=0; i<=j; i++)
		{
#ifdef BLAS
			clv[j][i] = cblas_ddot(j-i+1, &projStore[i][0], 1, &clv[j][i], 1);
#else
			reell Cnext = 0;
			for (int k=i; k<=j; k++)
				Cnext += projStore[i][k-i]*clv[j][k];
			clv[j][i] = Cnext;
#endif
		}



	//! 3) normalize coefficient matrix
	normalizeCLV(Ncol);

}




