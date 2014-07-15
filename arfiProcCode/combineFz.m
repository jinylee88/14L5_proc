function data = combineFz(lat, axial, arfidata, par, crossfadeWidthMm);
% function data = combineFz(lat, axial, arfidata, par, crossfadeWidthMm);
%
% Combine multiple focal zone data by averaging them together for each timestep
% 
% INPUTS:
%   lat - lateral axis (mm)
%   arfidata - ARFI displacement data (um) (designed to work with 3D data - axial x lat x t)
%   par - parameters structure from par file
%   crossfadeWidthMm - distance over which to blend the focal zones
%   
% OUTPUTS:
%   data - ARFI displacement data w/ focal zones merged together
%
% 
% Stephen Rosenzweig (stephen.rosenzweig@duke.edu)
% 2014-01-24

% pull out push focal depths
foci = par.pushFocalDepth;

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

% figure out what depth to blend the zones at
% average over lateral lines
meanData = squeeze(mean(arfidata,2));
% low-pass filter data
lambda = (par.c * 1e-3) / par.pushFreq;
[B A] = butter(4, double(1./(6*lambda./mean(diff(axial)))));
B = double(B);A = double(A);
meanData = filtfilt(B,A,double(meanData));
% crop the depth range to use to shallow to first focal depth and 
axInd = axial>(foci(1)-8*par.pushFnum^2*lambda) & axial<(foci(end)+8*par.pushFnum^2*lambda);

% tmpInd = 0; % debugging code
% average data within the crossfade regions
data = nan(size(arfidata,1), size(arfidata,2), size(arfidata,4));
for tInd = 1:size(arfidata,4);
    for i = 1:length(foci)
        if i==1
            axInd = axial>(foci(1)-8*par.pushFnum^2*lambda) & axial<(foci(2)+2*par.pushFnum^2*lambda);
            k = find(meanData(axInd,2,tInd)>=meanData(axInd,1,tInd), 1, 'first');
            if k>1
                crossfadeDepthMm = [axial(1), axial(k+find(axInd==1,1))-crossfadeWidthMm/2, axial(k+find(axInd==1,1))+crossfadeWidthMm/2];
            else
                % if you cannot find a good point where they cross, use the mean of the focal depths...probably need a better algorithm than this in
                % the future
                crossfadeDepthMm = [axial(1), mean(foci(i:i+1))-crossfadeWidthMm/2, mean(foci(i:i+1))+crossfadeWidthMm/2];
            end
        elseif i==length(foci)
            axInd = axial>(foci(end-1)-8*par.pushFnum^2*lambda) & axial<(foci(end)+2*par.pushFnum^2*lambda);
            k = find(meanData(axInd,end,tInd)>=meanData(axInd,end-1,tInd), 1, 'first');
            if k>1
                crossfadeDepthMm = [axial(k+find(axInd==1,1))+crossfadeWidthMm/2, axial(end)];
            else
                % if you cannot find a good point where they cross, use the mean of the focal depths...probably need a better algorithm than this in
                % the future
                crossfadeDepthMm = [mean(foci(i-1:i))+crossfadeWidthMm/2, axial(end)];
            end
        else
            axInd1 = axial>(foci(i-1)-8*par.pushFnum^2*lambda) & axial<(foci(i)+2*par.pushFnum^2*lambda);
            k1 = max(find(meanData(axInd1,i,tInd)>=meanData(axInd1,1,tInd), 1, 'first'),1);
            axInd2 = axial>(foci(i)-8*par.pushFnum^2*lambda) & axial<(foci(i+1)+2*par.pushFnum^2*lambda);
            k2 = max(find(meanData(axInd2,i+1,tInd)>=meanData(axInd2,i,tInd), 1, 'first'),1);
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

        data(crossfadeInd(1):crossfadeInd(2)-1,:,tInd) = arfidata(crossfadeInd(1):crossfadeInd(2)-1,:,i,tInd);
        % create linear profiles
        for j = 2:length(crossfadeInd)-1
            profileWeights = linspace(0,1,crossfadeInd(j+1)-crossfadeInd(j));
            profileWeights = profileWeights(:);
            tmp = squeeze(arfidata(crossfadeInd(j):crossfadeInd(j+1)-1,:,i,tInd)) .* repmat(profileWeights(end:-1:1), [1,size(data,2)]) + ...
                squeeze(arfidata(crossfadeInd(j):crossfadeInd(j+1)-1,:,i+1,tInd)) .* repmat(profileWeights, [1,size(data,2)]);
            data(crossfadeInd(j):crossfadeInd(j+1)-1,:,tInd) = tmp;
        end
        if i == length(foci)
            data(crossfadeInd(end-1):crossfadeInd(end),:,tInd) = arfidata(crossfadeInd(end-1):crossfadeInd(end),:,i,tInd);
        end

        clear crossfadeDepthMm crossfadeInd axInd axInd1 axInd2 k k1 k2
    end
end

