function Data = readDataOut(fileName)

    disp(['reading result netcdf file: ' fileName])

    %%%%%%%%%%%%%%%%%%%%%%%
    %% IMPORT OF NETcdf
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Get file ID:
    ncid = netcdf.open(fileName,'NC_NOWRITE');

    % Get variable IDs
    varid01= netcdf.inqVarID(ncid,'SW');
    varid02= netcdf.inqVarID(ncid,'SC');
    varid03= netcdf.inqVarID(ncid,'TW');
    varid04= netcdf.inqVarID(ncid,'TC');
    varid05= netcdf.inqVarID(ncid,'rateW');
    varid06= netcdf.inqVarID(ncid,'rateC');
    
    % Get the value of the variables, given their IDs.
    Data.SW = netcdf.getVar(ncid,varid01);
    Data.SC = netcdf.getVar(ncid,varid02);
    Data.TW = netcdf.getVar(ncid,varid03);
    Data.TC = netcdf.getVar(ncid,varid04);
    Data.rateW = netcdf.getVar(ncid,varid05);
    Data.rateC = netcdf.getVar(ncid,varid06);
    
    
    %get the dimensions to check whether this was calculated
    
    
    
    %read spike train if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'finalStatesSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no final states and currents stored in the netcdf file');
    end

    if (len > 0)
        disp('reading final states and currents');
        
        varidCur = netcdf.inqVarID(ncid, 'finalCurrents');
        Data.finalCurrents = netcdf.getVar(ncid, varidCur);
        
        varidFinal= netcdf.inqVarID(ncid, 'finalStates');
        finalState = netcdf.getVar(ncid, varidFinal);
        
        varidFinalIdx= netcdf.inqVarID(ncid, 'finalStatesIdx');
        finalStateIdx = netcdf.getVar(ncid, varidFinalIdx);
        
        %reshape the final state
        Data.finalStates = cell(1, length(finalStateIdx));
        finalStateIdx = finalStateIdx + 1; %c++ to matlab index conversion
        
        for f = 1:(length(finalStateIdx)-1)

            Data.finalStates(f) = {finalState((finalStateIdx(f)):(finalStateIdx(f+1)-1))};
            
        end
        Data.finalStates(end) = {finalState((finalStateIdx(end)):end)};
        
    end
    
    
    
    
    
    %read spike train if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'spikeTrainSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no spike train stored in the netcdf file');
    end

    if (len > 0)
        disp('reading spike train');
        varid07= netcdf.inqVarID(ncid, 'trainTime');
        varid08= netcdf.inqVarID(ncid, 'trainNeuron');
        Data.trainTime = netcdf.getVar(ncid, varid07);
        Data.trainNeuron = netcdf.getVar(ncid, varid08);
    end
    
    
    
    
    
    
    %read Lyapunov exponents if defined
    lenLE = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'LEonsSz');
        [name, lenLE] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no Lyapunov exponents stored in the netcdf file');
    end
    
    if (lenLE > 0)
        disp('reading Lyapunov exponents');
        varid09= netcdf.inqVarID(ncid, 'LEons');
        Data.LyapunovExponents = netcdf.getVar(ncid, varid09);
    
       
        
        
        %read times at which the orthonormalization were done
        lenTimes = 0;
        try
            dimid = netcdf.inqDimID(ncid, 'LEtimesSz');
            [name, lenTimes] = netcdf.inqDim(ncid, dimid);
        catch
            disp('no Lyapunov exponent times stored in the netcdf file');
        end

        if (lenTimes > 0)
            disp('reading Lyapunov exponents times');
            varid09b= netcdf.inqVarID(ncid, 'LEtimes');
            Data.LEtimes = netcdf.getVar(ncid, varid09b);
            
        end
        
        %read Lyapunov exponents convergence if defined
        lenConv = 0;
        try
            dimid = netcdf.inqDimID(ncid, 'LEconvergenceSz');
            [name, lenConv] = netcdf.inqDim(ncid, dimid);
        catch
            disp('no Lyapunov exponents convergence stored in the netcdf file');
        end

        if (lenConv > 0)
            disp('reading Lyapunov exponents convergence');
            varid09c= netcdf.inqVarID(ncid, 'LEconvergence');
            Data.LEconvergence = netcdf.getVar(ncid, varid09c);
            
            Data.LEconvergence = reshape(Data.LEconvergence, lenLE, lenConv/lenLE)';
        end
        
        
        %read backward iteration Lyapunov exponents if defined
        len = 0;
        try
            dimid = netcdf.inqDimID(ncid, 'LEclvSz');
            [name, len] = netcdf.inqDim(ncid, dimid);
        catch
            disp('no backward iteration Lyapunov exponents stored in the netcdf file');
        end

        if (len > 0)
            disp('reading backward iteration Lyapunov exponents');
            varid09d = netcdf.inqVarID(ncid, 'LEclv');
            Data.LEclv = netcdf.getVar(ncid, varid09d);
        end
        
        %read local Lyapunov exponents if defined
        len = 0;
        try
            dimid = netcdf.inqDimID(ncid, 'localLESz');
            [name, len] = netcdf.inqDim(ncid, dimid);
        catch
            disp('no local Lyapunov exponents stored in the netcdf file');
        end

        if (len > 0)
            disp('reading local Lyapunov exponents');
            varid09d = netcdf.inqVarID(ncid, 'localLE');
            Data.localLE = netcdf.getVar(ncid, varid09d);
            
            Data.localLE = reshape(Data.localLE, lenLE, len/lenLE)';
        end
    end
    
    
    
    
    
    
    %read firing rates if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'rateNeuronsSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no firing rates stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading firing rates');
        varid10= netcdf.inqVarID(ncid, 'rateNeurons');
        Data.rateNeurons = netcdf.getVar(ncid, varid10);
    end
    
    %read rate distribution
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'rateDistSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no rate distribution stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading rate distribution');
        varid11= netcdf.inqVarID(ncid, 'rateDistX');
        varid12= netcdf.inqVarID(ncid, 'rateDistY');
        Data.rateDistX = netcdf.getVar(ncid, varid11);
        Data.rateDistY = netcdf.getVar(ncid, varid12);
    end
    
    
    
    
    
    %read coefficients of variation if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'cvNeuronsSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no coefficient of variation stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading coefficients of variation');
        varid10= netcdf.inqVarID(ncid, 'cvNeurons');
        Data.cvNeurons = netcdf.getVar(ncid, varid10);
    end
    
    %read cv distribution
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'cvDistSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no coefficient of variation distribution stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading coefficient of variation distribution');
        varid11= netcdf.inqVarID(ncid, 'cvDistX');
        varid12= netcdf.inqVarID(ncid, 'cvDistY');
        Data.cvDistX = netcdf.getVar(ncid, varid11);
        Data.cvDistY = netcdf.getVar(ncid, varid12);
    end
    
    
    
    
        
    
    %read skewness if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'skewnessNeuronsSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no skewness stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading skewness');
        varid10= netcdf.inqVarID(ncid, 'skewnessNeurons');
        Data.skewnessNeurons = netcdf.getVar(ncid, varid10);
    end
    
    %read skewness distribution
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'skewnessDistSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no skewness distribution stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading skewness distribution');
        varid11= netcdf.inqVarID(ncid, 'skewnessDistX');
        varid12= netcdf.inqVarID(ncid, 'skewnessDistY');
        Data.skewnessDistX = netcdf.getVar(ncid, varid11);
        Data.skewnessDistY = netcdf.getVar(ncid, varid12);
    end
    
    
    
    
    
    
    
        
    
    %read kurtosis if defined
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'kurtosisNeuronsSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no kurtosis stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading kurtosis');
        varid10= netcdf.inqVarID(ncid, 'kurtosisNeurons');
        Data.kurtosisNeurons = netcdf.getVar(ncid, varid10);
    end
    
    %read cv distribution
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'kurtosisDistSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no kurtosis distribution stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading kurtosis distribution');
        varid11= netcdf.inqVarID(ncid, 'kurtosisDistX');
        varid12= netcdf.inqVarID(ncid, 'kurtosisDistY');
        Data.kurtosisDistX = netcdf.getVar(ncid, varid11);
        Data.kurtosisDistY = netcdf.getVar(ncid, varid12);
    end
    
    
    %read phase variables
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'phaseNeuronsSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no phases stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading phases');
        
        varid13= netcdf.inqVarID(ncid, 'phaseTimes');
        Data.phaseTimes = netcdf.getVar(ncid, varid13);
        
        varid14= netcdf.inqVarID(ncid, 'phaseNeurons');
        Data.phaseNeurons = netcdf.getVar(ncid, varid14);
        
        Data.phaseNeurons = reshape(Data.phaseNeurons, length(Data.phaseNeurons)/length(Data.phaseTimes), length(Data.phaseTimes));
    end
    
    %read distance
    len = 0;
    try
        dimid = netcdf.inqDimID(ncid, 'distancesSz');
        [name, len] = netcdf.inqDim(ncid, dimid);
    catch
        disp('no distances stored in the netcdf file');
    end
    
    if (len > 0)
        disp('reading distances');
        
        varid13= netcdf.inqVarID(ncid, 'phaseTimes');
        Data.phaseTimes = netcdf.getVar(ncid, varid13);
        
        varid14= netcdf.inqVarID(ncid, 'distances');
        Data.distances = netcdf.getVar(ncid, varid14);
        
    end
    
    
    netcdf.close(ncid);
    
end