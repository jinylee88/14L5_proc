function data = applyDdg_combineFz(lat, axial, arfidata, par, crossfadeWidthMm, normProfile);
% function [arfi]=combineFz(lat,arfidata,par);
%
% this function combines the separate applyDdg.m and combineFz.m functions with the change that the normalization profiles are used to determine the
% location of the merging of the focal depths
%
% the order of operations is: determine merging depths, apply DDG, merge data
% 
% INPUTS:
%   lat - lateral axis (mm)
%   axial - depth axis (mm)
%   arfidata - ARFI displacement data (um) (designed to work with 3D data - axial x lat x t)
%   par - parameters structure from par file
%   crossfadeWidthMm - distance over which to blend the focal zones
%   normProfile - normalization profile (axial x focal zone x time) (optional input)
%   
% OUTPUTS:
%   data - ARFI displacement data w/ focal zones merged together
%
% 
% Stephen Rosenzweig (stephen.rosenzweig@duke.edu)
% 2014-01-24

% pull out push focal depths
foci = par.pushFocalDepth;

% specify wavelength of push
lambda = (par.c * 1e-3) / par.pushFreq;

% baseline check for using separate focal zones
if ~par.separateFocalZoneAcqs
    error('Separate focal zone acquisitions not performed')
end

% check for crossfadeWidthMm distance making sense
if crossfadeWidthMm>min(diff(foci))/2
    warning('Crossfade width is greater than half the distance between pushes, ignoring regions where >2 foci overlap')
end

% make sure that there are the correct number of lateral positions for the number of focal zones specified
if((length(lat)*length(foci)) ~= size(arfidata,2)),
    error('%d lateral vectors present - %d lateral vectors expected', size(arfidata,2),(length(lat)*length(foci)));
end;

% reshape arfidata to isolate focal zones
arfidata = reshape(arfidata, size(arfidata,1), size(arfidata,2)/length(foci), length(foci), size(arfidata,3));

% if normalization profile is not specified, use average across field of view
if ~exist('normProfile', 'var') || isempty(normProfile)
    % average over lateral lines
    normProfile = squeeze(mean(arfidata,2));
    % low-pass filter data
    [B A] = butter(4, double(1./(6*lambda./mean(diff(axial)))));
    B = double(B);A = double(A);
    normProfile = filtfilt(B,A,double(normProfile));
end

% tmpInd = 0; % debugging code
% average data within the crossfade regions
data = nan(size(arfidata,1), size(arfidata,2), size(arfidata,4));
for tInd = 1:size(arfidata,4);
    for i = 1:length(foci)
        if i==1
            axInd = axial>(foci(1)-8*par.pushFnum^2*lambda) & axial<(foci(2)+2*par.pushFnum^2*lambda);
            k = find(normProfile(axInd,2,tInd)>=normProfile(axInd,1,tInd), 1, 'first');
            if k>1
                crossfadeDepthMm = [axial(1), axial(k+find(axInd==1,1))-crossfadeWidthMm/2, axial(k+find(axInd==1,1))+crossfadeWidthMm/2];
            else
                % if you cannot find a good point where they cross, use the mean of the focal depths...probably need a better algorithm than this in
                % the future
                crossfadeDepthMm = [axial(1), mean(foci(i:i+1))-crossfadeWidthMm/2, mean(foci(i:i+1))+crossfadeWidthMm/2];
            end
            normData(:,:,i) = arfidata(:,:,i,tInd)./repmat(normProfile(:,i,tInd),[1 size(arfidata,2)]);
            normData(:,:,i+1) = arfidata(:,:,i+1,tInd)./repmat(normProfile(:,i+1,tInd),[1 size(arfidata,2)]);
        elseif i==length(foci)
            axInd = axial>(foci(end-1)-8*par.pushFnum^2*lambda) & axial<(foci(end)+2*par.pushFnum^2*lambda);
            k = find(normProfile(axInd,end,tInd)>=normProfile(axInd,end-1,tInd), 1, 'first');
            if k>1
                crossfadeDepthMm = [axial(k+find(axInd==1,1))+crossfadeWidthMm/2, axial(end)];
            else
                % if you cannot find a good point where they cross, use the mean of the focal depths...probably need a better algorithm than this in
                % the future
                crossfadeDepthMm = [mean(foci(i-1:i))+crossfadeWidthMm/2, axial(end)];
            end
            normData(:,:,i-1) = arfidata(:,:,i-1,tInd)./repmat(normProfile(:,i-1,tInd),[1 size(arfidata,2)]);
            normData(:,:,i) = arfidata(:,:,i,tInd)./repmat(normProfile(:,i,tInd),[1 size(arfidata,2)]);
        else
            axInd1 = axial>(foci(i-1)-8*par.pushFnum^2*lambda) & axial<(foci(i)+2*par.pushFnum^2*lambda);
            k1 = max(find(normProfile(axInd1,i,tInd)>=normProfile(axInd1,1,tInd), 1, 'first'),1);
            axInd2 = axial>(foci(i)-8*par.pushFnum^2*lambda) & axial<(foci(i+1)+2*par.pushFnum^2*lambda);
            k2 = max(find(normProfile(axInd2,i+1,tInd)>=normProfile(axInd2,i,tInd), 1, 'first'),1);
            if isempty(k1),k1 = 1;end
            if isempty(k2),k2 = 1;end
            if k1>1 && k2>1
                crossfadeDepthMm = [axial(k1+find(axInd1==1,1))+crossfadeWidthMm/2, axial(k2+find(axInd2==1,1))-crossfadeWidthMm/2, axial(k2+find(axInd2==1,1))+crossfadeWidthMm/2];
            elseif k1>1
                crossfadeDepthMm = [axial(k1+find(axInd1==1,1))+crossfadeWidthMm/2, mean(foci(i:i+1))-crossfadeWidthMm/2, mean(foci(i:i+1))+crossfadeWidthMm/2];
            elseif k2>1
                crossfadeDepthMm = [mean(foci(i-1:i))+crossfadeWidthMm/2, axial(k2+find(axInd2==1,1))-crossfadeWidthMm/2, axial(k2+find(axInd2==1,1))+crossfadeWidthMm/2];
            else
                % if you cannot find a good point where they cross, use the mean of the focal depths...probably need a better algorithm than this in
                % the future
                crossfadeDepthMm = [mean(foci(i-1:i))+crossfadeWidthMm/2, mean(foci(i:i+1))-crossfadeWidthMm/2, mean(foci(i:i+1))+crossfadeWidthMm/2];
            end
            normData(:,:,i-1) = arfidata(:,:,i-1,tInd)./repmat(normProfile(:,i-1,tInd),[1 size(arfidata,2)]);
            normData(:,:,i) = arfidata(:,:,i,tInd)./repmat(normProfile(:,i,tInd),[1 size(arfidata,2)]);
            normData(:,:,i+1) = arfidata(:,:,i+1,tInd)./repmat(normProfile(:,i+1,tInd),[1 size(arfidata,2)]);
        end
        for j = 1:length(crossfadeDepthMm)
            crossfadeInd(j) = find(axial>=crossfadeDepthMm(j),1,'first');
        end
        
        
        
%         debugging code
%         if tmpInd~=tInd
%             fprintf(1, 'tInd = %d\n', tInd);
%         end
%         fprintf(1, '\t%0.2f  ', crossfadeDepthMm);
%         fprintf(1, '\n');
%         tmpInd = tInd;

        data(crossfadeInd(1):crossfadeInd(2)-1,:,tInd) = normData(crossfadeInd(1):crossfadeInd(2)-1,:,i);
        % create linear profiles
        for j = 2:length(crossfadeInd)-1
            profileWeights = linspace(0,1,crossfadeInd(j+1)-crossfadeInd(j));
            profileWeights = profileWeights(:);
            tmp = squeeze(normData(crossfadeInd(j):crossfadeInd(j+1)-1,:,i)) .* repmat(profileWeights(end:-1:1), [1,size(data,2)]) + ...
                squeeze(normData(crossfadeInd(j):crossfadeInd(j+1)-1,:,i+1)) .* repmat(profileWeights, [1,size(data,2)]);
            data(crossfadeInd(j):crossfadeInd(j+1)-1,:,tInd) = tmp;
        end
        if i == length(foci)
            data(crossfadeInd(end-1):crossfadeInd(end),:,tInd) = normData(crossfadeInd(end-1):crossfadeInd(end),:,i);
        end

        clear crossfadeDepthMm crossfadeInd axInd axInd1 axInd2 normData k k1 k2
    end
end

