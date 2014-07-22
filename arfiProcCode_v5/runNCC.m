function [u cc Iup Qup] = runNCC(I, Q, interpFactor, kernelLength,searchRegion, axial, par)
% function [u Iup Qup] = runLoupas(I, Q, interpFactor, axial, par)
%
% Inputs: I - in-phase data
%         Q - quadrature data
%         interpFactor - upsampling factor
%         axial - axial vector (used for demodulation)
%         par - parameters structure generated by arfi_image
%
% Outputs: u - displacement matrix
%          Iup - upsampled in-phase data
%          Qup - upsampled quadrature data

% Setup parameters
fs = par.fs*1e6*interpFactor;
fc = par.trackParams.reqList.txPulseReqs.TxPulse.fm_mhz*1e6; % Hz
c = 1540; % m/s
kasai_scale = c/(2*pi);

D = size(I); 
D(1) = D(1).*interpFactor;
if interpFactor>1
    [Iup,Qup] = computeUpsampledIQdata(I,Q,interpFactor);
else
    Iup = I;
    Qup = Q;
end
Iup = reshape(Iup, D);
Qup = reshape(Qup, D);

bwHz = 1e6*par.F0MHz*(par.trackParams.reqList.txPulseReqs.TxPulse.target_bw_pct/100);
nyquist = 2*(fc + 0.5*bwHz);
if fs <= nyquist
    error('Upsampled Fs (%0.2f MHz) is less than Nyquist (%0.2f MHz).',fs*1e-6,nyquist*1e-6)
end



