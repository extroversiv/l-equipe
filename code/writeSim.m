function [Hash, fileName] = writeSim(Para, directory)


%% Please simulation parameters
if (~isfield(Para, 'SW') && ~isfield(Para, 'TW'))
    Para.SW = 10; % SWC min. number of spikes per neuron during warmup before calculation(SWC*N),
    Para.TW = 0;   %minimum simulation time in s
end

if (~isfield(Para, 'SC') && ~isfield(Para, 'TC'))
    Para.SC = 10; %-SC min. number of spikes per neuron during calculation (SC*N)
    Para.TC = 0;   %provide either SC or TC and either SW or TW
end

if (isfield(Para, 'SW') && ~isfield(Para, 'TW'))
    Para.TW = 0;
end

if (~isfield(Para, 'SW') && isfield(Para, 'TW'))
    Para.SW = 0;
end

if (isfield(Para, 'SC') && ~isfield(Para, 'TC'))
    Para.TC = 0;
end

if (~isfield(Para, 'SC') && isfield(Para, 'TC'))
    Para.SC = 0;
end




if ~isfield(Para, 'ONstep')
    Para.ONstep = 1;
end
if ~isfield(Para, 'SR')
    Para.SR = 10; % min. number of spikes per neuron during warmup in current search(SWR*N)
end
if ~isfield(Para, 'TR')
    Para.TR = 0;
end
if ~isfield(Para, 'rateWnt')
    Para.rateWnt = 0;     %set to 0, if the external currents jsut set the final rate
end
if ~isfield(Para, 'pR')
    Para.pR = 0.01;        %precision of the desired rate
end

if ~isfield(Para, 'seedONS')
    Para.seedONS = 1;
end

if ~isfield(Para, 'train')
    Para.train = 0;
end
if ~isfield(Para, 'ISIneurons')
    Para.ISIneurons = 0;
end
if ~isfield(Para, 'ISIstats')
    Para.ISIstats = 0;
end
if ~isfield(Para, 'ISIbins')
    Para.ISIbins = 0;
end    
if ~isfield(Para, 'LyapunovExp')
    Para.LyapunovExp = 0;
end
if ~isfield(Para, 'pLE')
    Para.pLE = 0;        %precision in the main calculation
end

if ~isfield(Para, 'LyapunovExpConvergence')
    Para.LyapunovExpConvergence = 0;
end
if ~isfield(Para, 'SWONS')
    Para.SWONS = 1;
end
if ~isfield(Para, 'addCur');
    Para.addCur = 0;
end
   
if ~isfield(Para, 'pertSize');
    Para.pertSize = 0;
end

if ~isfield(Para, 'pertSpike');
    Para.pertSpike = 0;
end

if ~isfield(Para, 'pertSynapse');
    Para.pertSynapse = 0;
end

if ~isfield(Para, 'distances');
    Para.distances = 0;
end

if ~isfield(Para, 'phases');
    Para.phases = 0;
end

if ~isfield(Para, 'phaseSpikes');
    Para.phaseSpikes = [];
end

if ~isfield(Para, 'phaseTimes');
    Para.phaseTimes = [];
end

if ~isfield(Para, 'saveFinalState');
    Para.saveFinalState = 0;
end

if ~isfield(Para, 'CLV');
    Para.CLV = 0;
end

if ~isfield(Para, 'SWCLV');
    Para.SWCLV = 0;
end




%current vector and time vector must be the same length
if (Para.addCur)    
    if (Para.addCurHomo)
        if (length(Para.addCurTime) ~= length(Para.addCurIext))
            error('current vector and time vector must be the same length');
        end
    else
        if (length(Para.addCurTime) ~= size(Para.addCurIext, 2)) || (length(Para.addCurNeurons) ~= size(Para.addCurIext, 1))
            error('the dimensions of the current array are not correct');
        end
    end
end

%instantaneous population firing rate bin size
if ~isfield(Para, 'instPopRateBinSize');
    Para.instPopRateBinSize = 0;
end


Hash = DataHash(Para);
fileName = [directory, 'ParaSimulation-', Hash, '.nc'];
% ParaSimulation goes in front to avoid ncdump error (bug)
% ncdump: name begins with space or control-character: 1

%% Open netCDF files.
if ~exist(fileName, 'file')
    
    disp(['writing simulation netcdf file: ' fileName])

    ncid = netcdf.create(fileName, 'share');
    
    %% Define the dimensions of the variables.
    dimid_1 = netcdf.defDim(ncid, 'one', 1);
    dimid_szISI = netcdf.defDim(ncid, 'szISI', length(Para.ISIneurons));
    dimid_szTrain = netcdf.defDim(ncid, 'szTrain', length(Para.train));
    
    if (Para.addCur)
        dimid_szCurNeurons = netcdf.defDim(ncid, 'szCurNeurons', length(Para.addCurNeurons));
        dimid_szCurTime = netcdf.defDim(ncid, 'szCurTime', length(Para.addCurTime));
        
        if (Para.addCurHomo)
            dimid_szCurIext = netcdf.defDim(ncid, 'szCurIext', length(Para.addCurIext));
        else
            % save the 2D array with the additional current for each
            % neuron, as a 1D array, where each neuron's current is
            % appended after each other
            dimid_szCurIext = netcdf.defDim(ncid, 'szCurIext', length(Para.addCurIext)*length(Para.addCurNeurons));
            Para.addCurIext = reshape(Para.addCurIext, 1, length(Para.addCurIext)*length(Para.addCurNeurons));
        end
    end
    
    if (Para.phases || Para.distances)
        dimid_szPhaseSpikes = netcdf.defDim(ncid, 'szPhaseSpikes', length(Para.phaseSpikes));
        dimid_szPhaseTimes = netcdf.defDim(ncid, 'szPhaseTimes', length(Para.phaseTimes));
    end
    
    if (Para.pertSize)
        dimid_szPertVector = netcdf.defDim(ncid, 'szPertVector', length(Para.pertVector));
    end
    
    VarSim02 = netcdf.defVar(ncid, 'train', 'int', dimid_szTrain);
    VarSim03 = netcdf.defVar(ncid, 'ISIstats', 'int', dimid_1);
    VarSim04 = netcdf.defVar(ncid, 'ISIneurons', 'int', dimid_szISI);
    VarSim05 = netcdf.defVar(ncid, 'ISIbins', 'int', dimid_1);
    VarSim06 = netcdf.defVar(ncid, 'LyapunovExp', 'int', dimid_1);
    VarSim07 = netcdf.defVar(ncid, 'seedONS', 'int', dimid_1);
    VarSim08 = netcdf.defVar(ncid, 'SR', 'double', dimid_1);
    VarSim09 = netcdf.defVar(ncid, 'SW', 'double', dimid_1);
    VarSim10 = netcdf.defVar(ncid, 'SC', 'double', dimid_1);
    VarSim11 = netcdf.defVar(ncid, 'TR', 'double', dimid_1);
    VarSim12 = netcdf.defVar(ncid, 'TW', 'double', dimid_1);
    VarSim13 = netcdf.defVar(ncid, 'TC', 'double', dimid_1);
    VarSim15 = netcdf.defVar(ncid, 'rateWnt', 'double', dimid_1);
    VarSim16 = netcdf.defVar(ncid, 'pR', 'double', dimid_1);
    VarSim17 = netcdf.defVar(ncid, 'ONstep', 'int', dimid_1);
    VarSim18 = netcdf.defVar(ncid, 'saveFinalState', 'int', dimid_1);
    VarSim19 = netcdf.defVar(ncid, 'LyapunovExpConvergence', 'int', dimid_1);
    VarSim20 = netcdf.defVar(ncid, 'SWONS', 'double', dimid_1);
    VarSim21 = netcdf.defVar(ncid, 'pLE', 'double', dimid_1);
    VarSim22 = netcdf.defVar(ncid, 'CLV', 'int', dimid_1);
    VarSim23 = netcdf.defVar(ncid, 'SWCLV', 'int', dimid_1);
    VarSim24 = netcdf.defVar(ncid, 'instPopRateBinSize', 'double', dimid_1);

    
    
    VarSim_Cur = netcdf.defVar(ncid, 'addCur', 'int', dimid_1);

    if (Para.addCur)
        VarSim_CurHomo = netcdf.defVar(ncid, 'addCurHomo', 'int', dimid_1);
        VarSim_CurNeurons = netcdf.defVar(ncid, 'addCurNeurons', 'int', dimid_szCurNeurons);
        VarSim_CurTime = netcdf.defVar(ncid, 'addCurTime', 'double', dimid_szCurTime);
        VarSim_CurIext = netcdf.defVar(ncid, 'addCurIext', 'double', dimid_szCurIext);
    end
    
    VarSim_pertSize = netcdf.defVar(ncid, 'pertSize', 'double', dimid_1);
    VarSim_pertSpike = netcdf.defVar(ncid, 'pertSpike', 'int', dimid_1);
    VarSim_pertSynapse = netcdf.defVar(ncid, 'pertSynapse', 'int', dimid_1);
    
    VarSim_dist = netcdf.defVar(ncid, 'distances', 'int', dimid_1);
    VarSim_phases = netcdf.defVar(ncid, 'phases', 'int', dimid_1);
    
    if (Para.phases || Para.distances)
        if (~isempty(Para.phaseSpikes))
            VarSim_phaseSpikes = netcdf.defVar(ncid, 'phaseSpikes', 'int', dimid_szPhaseSpikes);
         else if (~isempty(Para.phaseTimes))
            VarSim_phaseTimes = netcdf.defVar(ncid, 'phaseTimes', 'double', dimid_szPhaseTimes);
             end
        end
    end
    if (Para.pertSize)
        VarSim_pertVector = netcdf.defVar(ncid, 'pertVector', 'double', dimid_szPertVector);
    end
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    %% Write data to variable.
    netcdf.putVar(ncid, VarSim02, Para.train - 1);  % C rrays start at 0 versus matlab arrays at 1
    netcdf.putVar(ncid, VarSim03, Para.ISIstats);
    netcdf.putVar(ncid, VarSim04, Para.ISIneurons - 1);
    netcdf.putVar(ncid, VarSim05, Para.ISIbins);
    netcdf.putVar(ncid, VarSim06, Para.LyapunovExp);
    netcdf.putVar(ncid, VarSim07, Para.seedONS);
    netcdf.putVar(ncid, VarSim08, Para.SR);
    netcdf.putVar(ncid, VarSim09, Para.SW);
    netcdf.putVar(ncid, VarSim10, Para.SC);
    netcdf.putVar(ncid, VarSim11, Para.TR*1000);
    netcdf.putVar(ncid, VarSim12, Para.TW*1000);
    netcdf.putVar(ncid, VarSim13, Para.TC*1000);
    netcdf.putVar(ncid, VarSim15, Para.rateWnt/1000);
    netcdf.putVar(ncid, VarSim16, Para.pR);
    netcdf.putVar(ncid, VarSim17, Para.ONstep);
    netcdf.putVar(ncid, VarSim18, Para.saveFinalState);
    netcdf.putVar(ncid, VarSim19, Para.LyapunovExpConvergence);
    netcdf.putVar(ncid, VarSim20, Para.SWONS);
    netcdf.putVar(ncid, VarSim21, Para.pLE);
    netcdf.putVar(ncid, VarSim22, Para.CLV);
    netcdf.putVar(ncid, VarSim23, Para.SWCLV);
    netcdf.putVar(ncid, VarSim24, Para.instPopRateBinSize);

    
    netcdf.putVar(ncid, VarSim_Cur, Para.addCur);
    
    if (Para.addCur)
        netcdf.putVar(ncid, VarSim_CurHomo, Para.addCurHomo);
        netcdf.putVar(ncid, VarSim_CurNeurons, Para.addCurNeurons - 1);
        netcdf.putVar(ncid, VarSim_CurTime, Para.addCurTime*1000);
        netcdf.putVar(ncid, VarSim_CurIext, Para.addCurIext);
    end
    
    netcdf.putVar(ncid, VarSim_pertSize, Para.pertSize);
    netcdf.putVar(ncid, VarSim_pertSpike, Para.pertSpike);
    netcdf.putVar(ncid, VarSim_pertSynapse, Para.pertSynapse);

    netcdf.putVar(ncid, VarSim_dist, Para.distances);
    netcdf.putVar(ncid, VarSim_phases, Para.phases);

    if (Para.phases || Para.distances)
        if (~isempty(Para.phaseSpikes))
            netcdf.putVar(ncid, VarSim_phaseSpikes, Para.phaseSpikes);
        else if (~isempty(Para.phaseTimes))
            netcdf.putVar(ncid, VarSim_phaseTimes, Para.phaseTimes*1000);
            end
        end
    end
    if (Para.pertSize)
        netcdf.putVar(ncid, VarSim_pertVector, Para.pertVector);
    end
        

    %% close file
    netcdf.close(ncid)
end
