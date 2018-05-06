function [varargout] = load_netcdf(adress_name,varargin)
% Function to load max 3-dimensional variables from netcdf. FK

nVarargs = length(varargin);
%fprintf('%s \n',[char(num2str(nVarargs)),' variable(s) from the file "',adress_name,'" to be loaded:'])

nc = netcdf.open(adress_name,'nowrite');
for k = 1:nVarargs;
    var_name = varargin{k};
    %fprintf('  %s \n',['-> ',var_name])
    varid = netcdf.inqVarID(nc,var_name); 
    
    var = double(netcdf.getVar(nc, varid));

    if size(var,2)~=1 && size(var,3)==1 && size(var,4)==1 % 2D 
        var = permute(var,[2 1]);    
    end
    if size(var,2)~=1 && size(var,3)~=1 && size(var,4)==1 % 3D 
        var = permute(var,[3 2 1]);    
    end    
    if size(var,2)~=1 && size(var,3)~=1 && size(var,4)~=1 % 4D 
        var = permute(var,[4 3 2 1]);    
    end   
    
    varargout{k} = var;
end
netcdf.close(nc)
