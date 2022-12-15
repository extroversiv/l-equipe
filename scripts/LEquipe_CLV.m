clear all
addpath('../code/')
set(0, 'defaultaxesfontsize', 20);

%%%%%%%%%%%%%% define the parameters of the network here %%%%%%%%%%%

neuronType = 1; %neuron type

N = 200;        %number of neurons
K = 50;         %number of synapses per neuron
J0 = -1;        %coupling strength
f = 5;          %network-averaged firing rate in Hz
tauM = 10;      %membrane time constant

rap = 10;        %AP onset rapidness in case of rapid theta neurons
tauS = tauM/2;  %synaptic time constant in case of cLIF or twoDlinear

%%%%%%%%%%%%%%%%%%%%%%%%%% end of input %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the given neuron parameters
ParaNet.N = N;
ParaNet.NeuronType = neuronType;
ParaNet.rapidness = rap;
ParaNet.tauM = tauM;

TwoDlinear.alpha = 1;
TwoDlinear.beta = 0;
TwoDlinear.gamma = 0;
TwoDlinear.delta = 1;
TwoDlinear.Cw = 0;
TwoDlinear.tauS = tauS;
ParaNet.twoDlinear = TwoDlinear;

%% set the random graph with K synapses per neuron on average
rand('twister', 1);
[ParaTopo.post ParaTopo.row_length] = random_graph(K, N);

%% set synapstic coupling strength (sqrt(K) scaling for the balanced state)
ParaTopo.J = J0/sqrt(K);

%% set the parameters of the simulation
ParaSim.rateWnt = f;        % this is the wanted firing rate
% the external currents that yield the wanted firing rate can be well
% approximated by the balance equation f = -I0/(J0*tauM)
% then with the balanced state scaling we end up with Iext = sqrt(K)*I0
ParaNet.Iext = -J0*f/1000*tauM*sqrt(K);

ParaSim.SW = 100;           % number of spikes per neuron during warmup

ParaSim.train = 1:N;        % neurons, whose spike times are saved

%Lyapunov exponent parameters
if ParaNet.NeuronType < 10
    ParaSim.LyapunovExp = ParaNet.N;
else
    ParaSim.LyapunovExp = 2*ParaNet.N;
end
ParaSim.SC = 15;            % avg. number of spikes per neuron in the calculation
ParaSim.SWONS = 10;         % warmup of the ONSE
ParaSim.ONstep = 1;         % orthonormalization step size

%covariant Lyapunov vectors
ParaSim.CLV = 1;            % calculate the CLVs
ParaSim.SWCLV = 10;         % warmup of the CLVs (must be < ParaSim.SC)

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

subplot(2,3,3)
plot(Data.trainTime, Data.trainNeuron, '.', 'markersize', 5);
xlabel('time');
ylabel('neurons')
xlim([0 max(max(Data.trainTime))])
title('spike train');

subplot(2,3,6)
plot(1/ParaSim.LyapunovExp:1/ParaSim.LyapunovExp:1, Data.LyapunovExponents)
hold all;
plot(1/ParaSim.LyapunovExp:1/ParaSim.LyapunovExp:1, Data.LEclv)
ylabel ('\lambda_i ( s ^{ -1})');
xlabel('i / N')
title('Lyapunov spectra')
legend(['forward '; 'backward'], 'Location', 'Northeast');


subplot(2,3,[1:2 4:5])
time = repmat(Data.LEtimes(1:size(Data.localLE, 1)), 1, size(Data.localLE, 2));
index = repmat(1/size(Data.localLE, 2):1/size(Data.localLE, 2):1, size(Data.localLE, 1), 1);
pcolor(time, index, Data.localLE);
shading flat;
zlim([min(min(Data.localLE)) max(max(Data.localLE))]);
xlabel('time (ms)');
ylabel('i / N')
title('local Lyapunov exponents per spike')
colormap(bluewhitered);
colorbar;

