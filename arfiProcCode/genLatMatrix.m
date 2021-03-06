function lat = genLatMatrix(par)
% function lat = genLatMatrix(par)
% 
% Inputs: par - parameters structure loaded from parameters.mat file, which is generated by arfi_image
% Outputs: lat - lateral position vector/matrix
%
% Use this function for overlapping receive regions (rx spacing ~= 1)

txBeams = par.trackParams.txXyzGridParams.P;
rxBeams = par.trackParams.rxXyzDelayGridParams.beamDeltaP;
lat = repmat(txBeams, [1 length(rxBeams)])'+repmat(rxBeams, [length(txBeams) 1])';
% txSpacing = nanmean(squeeze(par.trackParams.txXyzGridParams.fociMm(1,:,:))./squeeze(par.trackParams.txXyzGridParams.P));
switch par.probeType
    case 'linear'
        txSpacing = nanmean(squeeze(par.trackParams.txXyzGridParams.fociMm(1,1,:))./squeeze(par.trackParams.txXyzGridParams.P)); % nanmean in case there is a beam located at 0
    case 'phased'
        % this is currently an approximation....needs fixing
        txSpacing = nanmean(par.trackParams.txXyzGridParams.thetaX_DEG./par.trackParams.txXyzGridParams.P);
    case 'curvilinear'
        txSpacing = 1;
end
lat = lat*txSpacing';

% convert to a vector if there is no overlap between beam groups
if (par.trackParams.rxMultibeamParams.beamPatternP(end)-par.trackParams.rxMultibeamParams.beamPatternP(1))<=1
    lat = lat(:);
end