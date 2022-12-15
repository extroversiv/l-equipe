function [Hash, fileName] = writeTopo(Para, directory)

if ~isfield(Para, 'HomogSynapses')
    Para.HomogSynapses = 1;           %0 = false, 1 = true
end


if ~isfield(Para, 'pSyn')
    Para.pSyn = 1;           %0 = false, 1 = true
end


Hash = DataHash(Para);
fileName = [directory, 'ParaTopology-', Hash, '.nc'];
% ParaTopology goes in front to avoid ncdump error (bug)
% ncdump: name begins with space or control-character: 1

%% Open netCDF files.
if ~exist(fileName, 'file')
    
    disp(['writing topology netcdf file: ' fileName])

    ncid = netcdf.create(fileName, 'share');
    
    %% Define the dimensions of the variables.
    dimid_1 = netcdf.defDim(ncid, 'one', 1);
    dimid_N = netcdf.defDim(ncid, 'N', length(Para.row_length));
    dimid_post = netcdf.defDim(ncid, 'max_elements', length(Para.post));
    
    
    %% Define new variables in the topology file.
    VarSyn_Homo = netcdf.defVar(ncid, 'HomogSynapses', 'int', dimid_1);
    VarTopology_post = netcdf.defVar(ncid, 'post', 'int', dimid_post);
    VarTopology_rowLength = netcdf.defVar(ncid, 'row_length', 'int', dimid_N);

    if (Para.HomogSynapses)
        VarTopology_J = netcdf.defVar(ncid, 'J', 'double', dimid_1);
        VarTopology_pSyn = netcdf.defVar(ncid, 'pSyn', 'double', dimid_1);
    else
        VarTopology_J = netcdf.defVar(ncid, 'J', 'double', dimid_post);
        VarTopology_pSyn = netcdf.defVar(ncid, 'pSyn', 'double', dimid_post);
    end
    
 
    netcdf.putAtt(ncid, VarTopology_post, 'row_length', 'row_length')
    if (~Para.HomogSynapses)
        netcdf.putAtt(ncid, VarTopology_J, 'row_length', 'row_length')
        netcdf.putAtt(ncid, VarTopology_pSyn, 'row_length', 'row_length')
    end
    
    %% Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    %% Write data to variable.
    netcdf.putVar(ncid, VarSyn_Homo, Para.HomogSynapses);
    netcdf.putVar(ncid,VarTopology_post, Para.post - 1);  %convert from matlab to C indices
    netcdf.putVar(ncid,VarTopology_J, Para.J);
    netcdf.putVar(ncid,VarTopology_pSyn, Para.pSyn);
    netcdf.putVar(ncid,VarTopology_rowLength, Para.row_length);

    %% close file
    netcdf.close(ncid)
    
end
