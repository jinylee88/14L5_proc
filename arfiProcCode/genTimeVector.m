function [t, txTypeIndex, pushPRF] = genTimeVector(par)
% function t = genTimeVector(par)
%
% generates time vector for use with the SC2000 ARFI processing code
% Takes into account variable PRI, long delay
% Sets t=0 to transmit event of first track after the last push
%
% inputs: par - parameters structure saved with arfi_image
% outpus: t - time vector (ms)
%         pushPRF - pulse repetition frequency of the push (Hz)
%

% start with t=0 at the end of the push and work backwards for references
if par.separateFocalZoneAcqs
    t0 = (1:par.npush).*par.priusec(2);
    txTypeIndex0 = -ones(size(t0));
    t0(end+1:end+par.nref) = t0(end) + (1:par.nref).*par.priusec(1);
    txTypeIndex0(end+1:end+par.nref) = zeros(1,par.nref);
else
    t0 = (1:par.npush*length(par.pushFocalDepth)).*par.priusec(2);
    txTypeIndex0 = -ones(size(t0));
    t0(end+1:end+par.nref) = t0(end) + (1:par.nref).*par.priusec(1);
    txTypeIndex0(end+1:end+par.nref) = zeros(1,par.nref);
end

% now set up tracking times
t1 = (0:par.ntrack(1)).*par.priusec(3);
for i = 2:length(par.ntrack)
    t1(end+1:length(t1)+par.ntrack(i)) = t1(end) + (1:par.ntrack(i)).*par.priusec(i+2);
end
txTypeIndex1 = ones(size(t1));

t = [-t0(end:-1:1) t1(1:end-1)]; % strip off the last value of t1 (this is offtime after the last track)
t = t*1e-3; % convert to ms
txTypeIndex = [txTypeIndex0(end:-1:1) txTypeIndex1(1:end-1)];
if size(t)~= size(txTypeIndex),error('t and txTypeIndex vectors have different size');end
pushPRF = 1000./(t(end)-t(1)+(t1(end)-t1(end-1)).*1e-3); % factor in the last value of t1 into the push PRF and convert to Hz