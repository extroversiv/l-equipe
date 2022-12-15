function [Hash, fileName] = writeNet(Para, directory)


%% set the default values if not defined in Para
if ~isfield(Para, 'N')
    Para.N = 100;
end

if ~isfield(Para, 'HomogNetwork')
    Para.HomogNetwork = 1;           %0 = false, 1 = true
end

if ~isfield(Para, 'NeuronType')
    Para.NeuronType = 1;
    if (~Para.HomogNetwork)
        Para.NeuronType = Para.NeuronType*ones(1, Para.N);
    end
end
%1 = rapid theta
%2 = LIF
%3 = type1type2
%10 = twoDlinear (RIF->cLif)
%11 = cLIF (gamma = 1/3)
%12 = cLIF (gamma = 1/2)
%13 = cLIF (gamma = 2)
%14 = cLIF (gamma = 3)
%
%(4 = QIF,)
%(5-9 = other phase neurons)


if ~isfield(Para, 'tauM')
    Para.tauM = 10;
    if (~Para.HomogNetwork)
        Para.tauM = Para.tauM*ones(1, Para.N);
    end
end

if ~isfield(Para, 'Iext')
    Para.Iext = 1;
    if (~Para.HomogNetwork)
        Para.Iext = Para.Iext*ones(1,Para.N);
    end
end

if ~isfield(Para, 'seedInit')
    Para.seedInit = 1;
end

if ~isfield(Para, 'init')
    %if a seed is given, generate the initial states of the neurons
    rand('twister', Para.seedInit);
    Para.init = cell(1, Para.N);
    switch (Para.NeuronType)
        case {1, 3}
            for n = 1:Para.N
                Para.init{n} = 2*pi*(rand() - 0.5);
            end

        case {2}
            for n = 1:Para.N
                Para.init{n} = rand();
            end

        case  {10, 11, 12, 13, 14}
            for n = 1:Para.N
                Para.init{n} = [rand() 0];
            end
    end
end


if ~isfield(Para, 'rapidness')
    Para.rapidness = 1; %only for rapid theta
    if (~Para.HomogNetwork)
        Para.rapidness = Para.rapidness*ones(1, Para.N);
    end
end

if ~isfield(Para, 'type1type2Para')
    Para.type1type2Para = 0; %only for rapid theta
    if (~Para.HomogNetwork)
        Para.type1type2Para = Para.type1type2Para*ones(1, Para.N);
    end
end

if ~isfield(Para, 'twoDlinear')
    TwoDlinear.alpha = 1;
    TwoDlinear.beta = 0;
    TwoDlinear.gamma = 0;
    TwoDlinear.delta = 1;
    TwoDlinear.Cw = 0;
    TwoDlinear.tauS = Para.tauM./2;

    Para.twoDlinear = TwoDlinear;
    
    if (~Para.HomogNetwork)
        Para.twoDlinear.alpha = Para.twoDlinear.alpha*ones(1, Para.N);
        Para.twoDlinear.beta = Para.twoDlinear.beta*ones(1, Para.N);
        Para.twoDlinear.gamma = Para.twoDlinear.gamma*ones(1, Para.N);
        Para.twoDlinear.delta = Para.twoDlinear.delta*ones(1, Para.N);
        Para.twoDlinear.Cw = Para.twoDlinear.Cw*ones(1, Para.N);
    end
end




Hash = DataHash(Para);
fileName = [directory, 'ParaNeurons-', Hash, '.nc'];
% ParaNeurons goes in front to avoid ncdump error (bug)
% ncdump: name begins with space or control-character: 1


%% Open netCDF files.
if ~exist(fileName, 'file')
    
    disp(['writing neuron netcdf file: ' fileName])

    ncid = netcdf.create(fileName, 'share');

    %% Define the dimensions of the variables.
    dimid_1 = netcdf.defDim(ncid, 'one', 1);
    dimid_N = netcdf.defDim(ncid, 'N', Para.N);
    
    %% Define new variables in the neuron file.
    VarNeuron_N = netcdf.defVar(ncid, 'N', 'int', dimid_1);
    VarNeuron_Homo = netcdf.defVar(ncid, 'HomogNetwork', 'int', dimid_1);

    if (~Para.HomogNetwork)
        VarNeuron_NeuronType = netcdf.defVar(ncid, 'NeuronType', 'int', dimid_N);
        VarNeuron_tauM = netcdf.defVar(ncid, 'tauM', 'double', dimid_N);
        VarNeuron_Iext = netcdf.defVar(ncid, 'Iext', 'double', dimid_N);
        
        VarNeuron_rap = netcdf.defVar(ncid, 'rapidness', 'double', dimid_N);
        
        VarNeuron_type1type2 = netcdf.defVar(ncid, 'type1type2Para', 'double', dimid_N);

        VarNeuron_a = netcdf.defVar(ncid, 'twoDlinear_alpha', 'double', dimid_N);
        VarNeuron_b = netcdf.defVar(ncid, 'twoDlinear_beta', 'double', dimid_N);
        VarNeuron_c = netcdf.defVar(ncid, 'twoDlinear_gamma', 'double', dimid_N);
        VarNeuron_d = netcdf.defVar(ncid, 'twoDlinear_delta', 'double', dimid_N);
        VarNeuron_Cw = netcdf.defVar(ncid, 'twoDlinear_Cw', 'double', dimid_N);
        VarNeuron_tauS = netcdf.defVar(ncid, 'twoDlinear_tauS', 'double', dimid_N);
    else
        VarNeuron_NeuronType = netcdf.defVar(ncid, 'NeuronType', 'int', dimid_1);
        VarNeuron_tauM = netcdf.defVar(ncid, 'tauM', 'double', dimid_1);
        VarNeuron_Iext = netcdf.defVar(ncid, 'Iext', 'double', dimid_1);
        
        VarNeuron_rap = netcdf.defVar(ncid, 'rapidness', 'double', dimid_1);
        
        VarNeuron_type1type2 = netcdf.defVar(ncid, 'type1type2Para', 'double', dimid_1);

        VarNeuron_a = netcdf.defVar(ncid, 'twoDlinear_alpha', 'double', dimid_1);
        VarNeuron_b = netcdf.defVar(ncid, 'twoDlinear_beta', 'double', dimid_1);
        VarNeuron_c = netcdf.defVar(ncid, 'twoDlinear_gamma', 'double', dimid_1);
        VarNeuron_d = netcdf.defVar(ncid, 'twoDlinear_delta', 'double', dimid_1);
        VarNeuron_Cw = netcdf.defVar(ncid, 'twoDlinear_Cw', 'double', dimid_1);
        VarNeuron_tauS = netcdf.defVar(ncid, 'twoDlinear_tauS', 'double', dimid_1);
    end

    
    state = [];
    stateIdx = zeros(1, length(Para.init));
    index = 0;
    for s = 1:length(Para.init)
        stateIdx(s) = index;
        state = [state Para.init{s}];
        index = index + length(Para.init{s});
    end
    
    dimid_state = netcdf.defDim(ncid, 'initStatesSz', length(state));
    dimid_stateIdx = netcdf.defDim(ncid, 'initStatesIdxSz', length(stateIdx));
    
    VarNeuron_state = netcdf.defVar(ncid, 'initStates', 'double', dimid_state);
    VarNeuron_stateIdx = netcdf.defVar(ncid, 'initStatesIdx', 'int', dimid_stateIdx);
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);


    %% Write data to variable.
    netcdf.putVar(ncid, VarNeuron_N, Para.N);
    netcdf.putVar(ncid, VarNeuron_Homo, Para.HomogNetwork);
    netcdf.putVar(ncid, VarNeuron_NeuronType, Para.NeuronType);
    netcdf.putVar(ncid, VarNeuron_tauM, Para.tauM);
    netcdf.putVar(ncid, VarNeuron_Iext, Para.Iext);
    netcdf.putVar(ncid, VarNeuron_state, state);
    netcdf.putVar(ncid, VarNeuron_stateIdx, stateIdx);
    
    netcdf.putVar(ncid, VarNeuron_rap, Para.rapidness);

    netcdf.putVar(ncid, VarNeuron_type1type2, Para.type1type2Para);

    netcdf.putVar(ncid, VarNeuron_a, Para.twoDlinear.alpha);
    netcdf.putVar(ncid, VarNeuron_b, Para.twoDlinear.beta);
    netcdf.putVar(ncid, VarNeuron_c, Para.twoDlinear.gamma);
    netcdf.putVar(ncid, VarNeuron_d, Para.twoDlinear.delta);
    netcdf.putVar(ncid, VarNeuron_Cw,   Para.twoDlinear.Cw);
    netcdf.putVar(ncid, VarNeuron_tauS, Para.twoDlinear.tauS);
   
    
    %% close file
    netcdf.close(ncid);
    
end