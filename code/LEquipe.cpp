/*
 * lequipe.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

/*
 * git add -A
 * git commit
 * git push
 *
 * Pull from repo
 * git remote add origin https://mig80@bitbucket.org/mig80/lequipe.git
 * git push -u origin master
 *
 * Undo a commit and redo
 * $ git commit ...
 * $ git reset --soft HEAD^      (1)
 * $ edit                        (2)
 * $ git add ....                (3)
 * $ git commit -c ORIG_HEAD     (4)
 * (1) This is most often done when you remembered what you just committed is incomplete, or you misspelled your commit message, or both. Leaves working tree as it was before "reset".
 * (2) Make corrections to working tree files.
 * (3) Stage changes for commit.
 * (4) "reset" copies the old head to .git/ORIG_HEAD; redo the commit by starting with its log message. If you do not need to edit the message further, you can give -C option instead.
 *
 * reset to commit cd50aeb71477 online:
 * git push origin +cd50aeb71477:master
 */

/*
 * git push
 */

#include <iostream>											//!< standard C++ input output library
#include <sstream>											//!< strings
#include <fstream>											//!< file handling
#include <limits>											//!< show limits of data types
#include <ctime>											//!< for time calculations
#include <cmath>											//!< for simple math functions
#include "reell.h"											//!< provides a simple change of double precision
#include "LEnetwork.h"
#include "LEanalysis.h"
#include "LEinputoutput.h"
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;										//!< standard C++ namespace



//! The main program. 
/*!
  The main C++ program that handles all simulations.
*/

int main(int argc, char** argv)
{

	int nP = 1;																	//!< number of processors involved
	int myID = 0;																//!< local ID of the processors, for unique communications

#ifdef PAR
	//! Initialize the parallel environment.

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nP);
	MPI_Comm_rank(MPI_COMM_WORLD,&myID);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//! Measure the CPU time
	time_t start, end, startAll;
	reell CPUtime;
	if (myID == 0)
	{
		time(&start);
		time(&startAll);
	}



	if (myID == 0) cout << endl << "running on " << nP << " processor(s) ..." << endl;


	st_in in;
	st_out out;

	LE_input_output inout;


#ifdef PAR
	//! the following is displayed successively by each processor
	int dummy = 0;
	MPI_Status status;

	if (myID > 0)
		MPI_Recv(&dummy, 1, MPI_INT, myID - 1, 0, MPI_COMM_WORLD, &status);			//!< receive dummy message from previous processor

#endif
	cout << endl << "************************ CPU " << myID + 1 << " of " << nP << " reporting for duty ************************" << endl;

	//! display the hostname and the current working directory
	char path[128];
	char hostname[128];

	gethostname(hostname, sizeof(hostname));
	getcwd(path, sizeof(path));

	cout << "running on " << hostname << " in " << path << endl;
	cout << "the " << argc - 1 << " command line arguments were: ";
	for (int i=1; i<argc; i++)
		cout << argv[i] << " ";
	cout << endl;

	cout << "-------------------------------------------------------------------------------" << endl;

	//! display the data types that will be used
	cout << "the used data types are:" << endl;
	cout << "int:\t\t size = " << sizeof(int) << " Bytes,\t range = " << numeric_limits<int>::min() << " to " << numeric_limits<int>::max() << endl;
	cout << "unsigned:\t size = " << sizeof(unsigned) << " Bytes,\t range = " << numeric_limits<unsigned>::min() << " to " << numeric_limits<unsigned>::max() << endl;
	cout << "long long:\t size = " << sizeof(long long) << " Bytes,\t range = " << numeric_limits<long long>::min() << " to " << numeric_limits<long long>::max() << endl;
	cout << "reell:\t\t size = " << sizeof(reell) << " Bytes,\t range = " << numeric_limits<reell>::min()
			<< " to " << numeric_limits<reell>::max() << ",\t precision = " << numeric_limits<reell>::epsilon() << endl;

	cout << "-------------------------------------------------------------------------------" << endl;

	if (argc < 5)
	{
		if (myID == 0) cout << " There need to be 4 filenames after the program name with the input and output filenames," << endl
				<< "like ./LEquipe 12_neurons.nc 34_topology.nc 56_simulation.nc 78_results.nc" << endl;
		return 1;
	}

	//! The input parameters in the command line are the 4 netcdf files, that must come in the following order
	//! 12_neurons.nc 34_topology.nc 56_simulation.nc 78_results.nc
	//! The last input parameter is the filename for the ouput file containing the results

	string fileNeurons = argv[1];
	string fileTopology = argv[2];
	string fileSimulation = argv[3];
	string fileResults = argv[4];

	inout.read_netcdf(fileNeurons, fileTopology, fileSimulation, &in);
	out.N = in.Nall;
	out.phases = in.phases;

	//! for the perturbed trajectory
	st_out outPert = out;


#ifdef PAR
	sleep(0.1);		//!< wait 0.1 second to synchronize the outputs from different processors

	if (myID < nP - 1)
		MPI_Send(&dummy, 1, MPI_INT, myID + 1, 0, MPI_COMM_WORLD);

	//! wait for all processor to finish
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	if (myID == 0)
	{
		time(&end);
		CPUtime = difftime(end, start);
		time(&start);
		cout << endl << "CPU time for initialization: " << CPUtime << "s" << endl;
	}

	if (myID == 0) cout << endl << "**************************** seting up the network ****************************" << endl;


	LE_network network(&in);
	LE_analysis simulation(&network);


	if (myID == 0)
	{
		time(&end);
		CPUtime = difftime(end, start);
		time(&start);
		cout << endl << "CPU time for network setup: " << CPUtime << "s" << endl;
	}

	if (myID == 0) cout << endl << "************************** starting the simulations ***************************" <<  endl;


	simulation.setRate_warmup(&in, &out);


	//! save the original state
	vector<vector<reell> > stateOrg = network.get_state_Nloc();

	bool applyPerturbation = false;
	simulation.multitask(&in, &out, applyPerturbation);

	if (myID == 0)
	{
		time(&end);
		CPUtime = difftime(end, start);
		time(&start);
		cout << endl << "CPU time for simulation: " << CPUtime << "s" << endl;


		cout << endl << "***************************** saving the results ******************************" << endl << endl;
	}


	inout.write_netcdf(fileResults, &out, fileNeurons, fileTopology, fileSimulation);




	if (in.pertSize || in.pertSpike || in.pertSynapse)
	{
		if (myID == 0) cout << endl << "***************************** perturbed trajectory ******************************" << endl << endl;

		//! reset to the original state
		network.set_state_Nloc(stateOrg);

		applyPerturbation = true;
		simulation.multitask(&in, &outPert, applyPerturbation);

		if (in.distances)
		{
			//! Calculate the distances with the norm provided in in.distance
			outPert.distances = vector<reell> (out.phaseTimes.size(), 0);
			for (unsigned t=0; t<out.phaseTimes.size(); t++)
				for (unsigned n=0; n<out.phaseNeurons[t].size(); n++)
					outPert.distances[t] += gsl_pow_int(abs(outPert.phaseNeurons[t][n] - out.phaseNeurons[t][n]), in.distances);

			int norm = out.phaseNeurons[0].size();

#ifdef PAR
			//! Gather sum from all processors on root
			if (myID == 0)
			{
				vector<double> distAll(out.phaseTimes.size());
				int normAll;

				MPI_Reduce(&out.distances.front(), &distAll.front(), out.phaseTimes.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				out.distances = distAll;

				MPI_Reduce(&norm, &normAll, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				norm = normAll;
			}
			else
			{
				double dummy;
				int dummy2;

				MPI_Reduce(&out.distances.front(), &dummy, out.phaseTimes.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&norm, &dummy2, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			}
#endif

			if (myID == 0)
				for (unsigned t=0; t<out.phaseTimes.size(); t++)
					outPert.distances[t] /= norm;
		}



		if (myID == 0)
		{
			time(&end);
			CPUtime = difftime(end, start);
			time(&start);
			cout << endl << "CPU time for simulation: " << CPUtime << "s" << endl;


			cout << endl << "***************************** saving the results ******************************" << endl << endl;
		}


		stringstream fileResultsPert;
		fileResultsPert << fileResults << "_perturbed";
		inout.write_netcdf(fileResultsPert.str(), &outPert, fileNeurons, fileTopology, fileSimulation);


	}


	if (myID == 0)
	{
		time(&end);
		CPUtime = difftime(end, start);
		time(&start);
		cout << endl << "CPU time for saving results: " << CPUtime << "s" << endl;

		cout << endl << "********************************** the end ************************************" << endl << endl;

		time(&end);
		CPUtime = difftime(end, startAll);
		cout << "CPU time overall: " << CPUtime << "s" << endl << endl;
	}

#ifdef PAR
	MPI_Finalize();
#endif

	return 0;
}
	
