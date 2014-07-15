function data = applyDdg(axial, arfidata, par, scaleDisp, normProfile);
% function data = applyDdg(axial, arfidata, normProfile);
%
% applys depth-dependent gain based on either an given normalization profile or the mean of the data across all lateral positions
%
% inputs: axial - axial position vector (mm)
%         arfidata - displacement data (um), (3D matrix - axial x lat x time)
%         par - parameters struct from par_*.mat
%         scaleDisp - whether (1) or not (0) to scale normalization profile to [0 1] (default is 0)
%         normProfile - normalization profile to use (optional, can be any of the following sizes)
%                       3D matrix - axial x focal zone x time
%                       2D matrix - axial x time
%                       2D matrix - axial x focal zone
%                       1D vector - axial
%                       default value is the mean across the field of view (of each focal depth)
%
% output: data - normalized arfidata

if ~exist('scaleDisp', 'var'),scaleDisp = 0;end

if ndims(arfidata)<2
    error('arfidata is incorrectly formed');
end

% low pass filter for normalization profiles
lambda = (par.c * 1e-3) / par.pushFreq;
[B A] = butter(4, double(1./(6*lambda./nanmean(diff(axial)))));
B = double(B);A = double(A);

if par.separateFocalZoneAcqs
    % pull out push focal depths
    foci = par.pushFocalDepth;
    % reshape arfidata to isolate focal zones
    arfidata = reshape(arfidata, size(arfidata,1), size(arfidata,2)/length(foci), length(foci), size(arfidata,3));
    % decide how to normalize data data depending on if normProfile was an input
    if ~exist('normProfile', 'var')
        normProfile = squeeze(nanmean(arfidata,2));
        normProfile = single(filtfilt(B,A,double(normProfile)));
        % put normProfile on a [0 1] scale
        if scaleDisp
            for i = 1:size(normProfile,2)
                axInd = axial<(foci(i)+4*par.pushFnum.^2*lambda); % don't look for max deep to the depth of field
                normProfile(:,i,:) = normProfile(:,i,:) ./ repmat(max(normProfile(axInd,i,:), [], 1), [size(normProfile,1), 1, 1]);
            end
        end
    elseif ndims(normProfile)==1
        normProfile = repmat(normProfile, [1, size(arfiata,3), size(arfidata,4)]);
    elseif ndims(normProfile)==2
        if size(normProfile,2) == size(arfidata,3) % axial x focal zone
            normProfile = repmat(normProfile, [1, 1, size(arfiata,4)]);
        elseif size(normProfile,2) == size(arfidata,4) % axial x time
            normProfile = permute(repmat(normProfile, [1, 1, size(arfiata,3)]), [1 3 2]);
        else
            error('Normalization profile is 2D, but does not match sizes of arfidata')
        end
    end
    % make sure that the sizes are correct
    if size(normProfile,1)~=size(arfidata,1) || size(normProfile,2)~=size(arfidata,3) || size(normProfile,3)~=size(arfidata,4)
        error('Normalization profile and arfidata dimension sizes do not match')
    end
    normProfile(normProfile<0.4) = 0.4; % cap the maximum amplification
    normProfile = single(filtfilt(B,A,double(normProfile))); % re-filter to remove sharp edges
    % apply normalization profile
    normProfile = permute(repmat(normProfile, [1,1,1,size(arfidata,2)]), [1 4 2 3]); % make normProfile the same size as arfidata
    data = arfidata./normProfile;
    data = reshape(data, size(data,1), [], size(data, 4));
else
    % decide how to normalize data data depending on if normProfile was an input
    if ~exist('normProfile', 'var')
        normProfile = squeeze(nanmean(arfidata,2));
        normProfile = single(filtfilt(B,A,double(normProfile)));
        % put normProfile on a [0 1] scale
        if scaleDisp
            axInd = axial<(max(par.pushFocalDepth)+4*par.pushFnum.^2*lambda); % don't look for max deep to the depth of field
            normProfile = normProfile ./ repmat(max(normProfile(axInd,:), [], 1), [size(normProfile,1), 1]);
        end
    elseif ndims(normProfile)==1
        normProfile = repmat(normProfile, [1 size(arfiata,3)]);
    end
    % make sure that the sizes are correct
    if size(normProfile,1)~=size(arfidata,1) || size(normProfile,2)~=size(arfidata,3)
        error('Normalization profile and arfidata dimension sizes do not match')
    end
    normProfile(normProfile<0.4) = 0.4; % cap the maximum amplification
    normProfile = single(filtfilt(B,A,double(normProfile))); % re-filter to remove sharp edges
    % apply normalization profile
    normProfile = permute(repmat(normProfile, [1 1 size(arfidata,2)]), [1 3 2]); % make normProfile the same size as arfidata
    data = arfidata./normProfile;
end