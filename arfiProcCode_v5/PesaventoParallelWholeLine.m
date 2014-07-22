function [Tau CC]= PesaventoParallelWholeLine(IQ0,IQ1,fs,f0,SrchWin,N)
if ~exist('SrchWin','var')
    SrchWin = 100e-6/770*fs;
end
if ~exist('N','var')
    N = 3;
end
grossLags = [-ceil(SrchWin):ceil(SrchWin)];
IQ0 = repmat(IQ0,size(IQ1)./size(IQ0));
tau = zeros([size(IQ0)]);
Tau = zeros(size(IQ0));
K = size(IQ0,1);
if nargout > 1
    CC = zeros(size(IQ0));
end

x0 = IQ0;
x1 = IQ1;

mux0 = nanmean(x0,1);
mux1 = nanmean(x1,1);
sigx0 = nanstd(x0,1,1);
sigx1 = nanstd(conj(x1),1,1);
x0bZeroMean = x0-mux0(ones(K,1),:,:,:);
x1bZeroMean = x1-mux1(ones(K,1),:,:,:);

sz = size(x0);
NCC = zeros([length(grossLags),sz(2:end)]);
for i = 1:length(grossLags);
    idx = max(1,min(K,(1:K)+grossLags(i)));
    NCC(i,:,:) = nansum(conj(x1bZeroMean(idx,:,:)).*x0bZeroMean)./((K).*(sigx0.*sigx1));
end
[CCpk CCpkidx] = max(NCC,[],1);
SampleShift = grossLags(CCpkidx);

sigx1SampleShifted = sigx1;
x1SampleShifted = x1;

for j = 1:size(x1,2)
    for k = 1:size(x1,3);
        idx = max(1,min(K,(1:K)+SampleShift(1,j,k)));
        x1SampleShifted(:,j,k) = x1(idx,j,k).*exp(1j*2*pi*f0*SampleShift(ones(K,1),j,k)*(1/fs));
    end
end
x1 = x1SampleShifted;

corrVal = nansum(conj(x1-repmat(nanmean(x1,1),[K 1 1])).*(x0-repmat(nanmean(x0,1),[K 1 1])))./((K).*sigx1.*sigx0);
corrValPhasCrct = angle(corrVal);
tau = corrValPhasCrct/(2*pi*f0);

for n = 1:N-1
    x0PhaseShifted = subsampleshift(x0,f0,tau/2);
    x1PhaseShifted = subsampleshift(x1,f0,-tau/2);
    corrVal = nansum(conj(x1PhaseShifted-repmat(nanmean(x1PhaseShifted,1),[K 1 1])).*(x0PhaseShifted-repmat(nanmean(x0PhaseShifted,1),[K 1 1])))./((K).*sigx1.*sigx0);
    corrValPhasCrct = angle(exp(-1j*2*pi*f0*(tau)/fs).*corrVal);
    tau = tau + corrValPhasCrct/(2*pi*f0); 
end
if nargout > 1
    CC = squeeze(corrVal);
end

Tau = tau+SampleShift*(1/fs);
end

