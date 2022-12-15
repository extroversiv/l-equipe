clear all
addpath('../code/')
set(0, 'defaultaxesfontsize', 20);

%%%%%%%%%%%%%% define the parameters of the network here %%%%%%%%%%%

neuronType = 1; %neuron type

NI = 200;        %number of inhibitory neurons
NE = 800;        %number of excitatory neurons
K = 50;         %number of synapses per neuron
J0 = 1;        %coupling strength
f = 5;          %network-averaged firing rate in Hz
tauM = 10;      %membrane time constant

rap = 10;        %AP onset rapidness in case of rapid theta neurons
tauS = tauM/2;  %synaptic time constant in case of cLIF or twoDlinear

%%%%%%%%%%%%%%%%%%%%%%%%%% end of input %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the given neuron parameters
N = NE + NI;

ParaNet.HomogNetwork = 0;
ParaNet.N = N;
ParaNet.NeuronType = neuronType*ones(N,1);
ParaNet.rapidness = rap*ones(N,1);
ParaNet.tauM = tauM*ones(N,1);

TwoDlinear.alpha = 1*ones(N,1);
TwoDlinear.beta = 0*ones(N,1);
TwoDlinear.gamma = 0*ones(N,1);
TwoDlinear.delta = 1*ones(N,1);
TwoDlinear.Cw = 0*ones(N,1);
TwoDlinear.tauS = tauS*ones(N,1);
ParaNet.twoDlinear = TwoDlinear;

%% create the mixed network with K synapses per neuron on averageParaNet.twoDlinear.Cw = ParaNet.twoDlinear.Cw*ones(N,1);

%% set synapstic coupling strength (sqrt(K) scaling for the balanced state)
rand('twister', 1);

A_EE = 0.01*J0/sqrt(K)*(rand(NE) < K/NE);
A_II = -J0/sqrt(K)*(rand(NI) < K/NI);

A_IE = 0.01*J0/sqrt(K)*(rand(NE, NI) < K/NI);
A_EI = -J0/sqrt(K)*(rand(NI, NE) < K/NE);

A = horzcat(vertcat(A_EE, A_EI), vertcat(A_IE, A_II));

[ParaTopo.post, ParaTopo.J, ParaTopo.row_length] = adjacencymatrix2postsynapticvector(A);

ParaTopo.HomogSynapses = 0;
ParaTopo.pSyn = ones(size(ParaTopo.J));

ParaNet.Iext = vertcat(ones(NE,1), 0.00001*ones(NI,1));
ParaNet.init = num2cell(2*pi*(rand(N,1) - 0.5));

%% set the parameters of the simulation
ParaSim.rateWnt = f;        % this is the wanted firing rate
% the external currents that yield the wanted firing rate can be well
% approximated by the balance equation f = -I0/(J0*tauM)
% then with the balanced state scaling we end up with Iext = sqrt(K)*I0

ParaSim.SW = 100;           % number of spikes per neuron during warmup
ParaSim.TC = 1;             % time duration of the calculation in seconds

ParaSim.train = 1:N;        % neurons, whose spike times are saved


ParaSim.LyapunovExp = N;
ParaSim.ONstep = 10*N/K;

%% write all parameters to netcdf files to directory data/ and get the hashes
directory = '../data/';
if ~exist(directory, 'dir')
    disp(['creating new directory: ' directory]);
    mkdir(directory)
end

[HashNet, FileNet] = writeNet(ParaNet, directory);
[HashTopo, FileTopo] = writeTopo(ParaTopo, directory);
[HashSim, FileSim] = writeSim(ParaSim, directory);
HashDataOut = DataHash([HashNet, HashTopo, HashSim]);
FileOut = [directory, 'DataOut-', HashDataOut, '.nc'];

%% run the C++ simulation
system(['../LEquipe ', FileNet, ' ', FileTopo, ' ', FileSim, ' ', FileOut]);

%% read the output file and plot the results
Data = readDataOut(FileOut);


figure;
subplot(2,2,1)
plot(Data.trainTime, Data.trainNeuron, '.', 'markersize', 5);
xlabel('time (s)');
ylabel('neurons')
xlim([0 ParaSim.TC])

subplot(2,2,2)
plot(1/ParaSim.LyapunovExp:1/ParaSim.LyapunovExp:1, Data.LyapunovExponents)
ylabel ('\lambda_i ( s ^{ -1})');
xlabel('i / N')

