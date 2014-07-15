function data = ccFiltData(lat, axial, arfidata, cc_coef, ccThreshold);
% function data = ccFiltData(lat, axial, arfidata, cc_coef, ccThreshold);
%
% remove data with cc_coef values less than ccThreshold
% then re-interpolate using TriScatteredInterp
%
% inputs: lat - lateral position vector (mm)
%         axial - axial position vector (mm)
%         arfidata - displacement data (um), (3D matrix - axial x lat x time)
%         cc_coef - correlation coefficients (0-1 scale), (3D matrix - axial x lat x time)
%         ccThreshold - correlation coefficient threshold (default value = 0.95)
%
% output: data - re-interpolated data, same size as arfidata

if ~exist('ccThreshold', 'var'),ccThreshold = 0.95;end % default value for ccThreshold

% make X,Z matrices
lat = (1:size(arfidata,2)).*mean(diff(lat));lat = lat-mean(lat); % doing this for the separateFocalZones cases where length(lat)~=size(arfidata,2)
[X,Z] = meshgrid(lat, axial);

data = nan(size(arfidata), 'single'); % initialize matrix
for ts = 1:size(arfidata,3)
    
    % output which timestep is currently being filtered
    if ts==1
        fprintf(1, 'Correlation coefficient filter for time step %d/%d', ts, size(arfidata,3));
    elseif ts==size(arfidata,3)
        tmpS = sprintf('%d/%d', ts-1, size(arfidata,3));
        fprintf(1, repmat('\b', [1 length(tmpS)]));
        fprintf(1, '%d/%d', ts, size(arfidata,3));
        fprintf(1, '\n');
    else
        tmpS = sprintf('%d/%d', ts-1, size(arfidata,3));
        fprintf(1, repmat('\b', [1 length(tmpS)]));
        fprintf(1, '%d/%d', ts, size(arfidata,3));
    end
    
    % perform filtering and interpolation operations
    tmp = arfidata(:,:,ts);
    mask = cc_coef(:,:,ts)>ccThreshold; % figure out which data to filter out
    X1 = X(mask);
    Z1 = Z(mask);
    tmp = tmp(mask);
    F = TriScatteredInterp(double(X1),double(Z1),double(tmp));
    out = F(X(:),Z(:));
    out = reshape(out, size(X));
    out(isnan(out)) = 0; % remove the nan values and set to 0 instead
    data(:,:,ts) = single(out);
end