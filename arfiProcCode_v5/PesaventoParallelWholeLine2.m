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

mux0 = nanmean(reshape(x0,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]),1);
mux1 = nanmean(reshape(x1,[size(x1,1)*size(x1,2),size(x1,3),size(x1,4)]),1);
sigx0 = nanstd(reshape(x0,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]),1,1);
sigx1 = nanstd(conj(reshape(x1,[size(x1,1)*size(x1,2),size(x1,3),size(x1,4)])),1,1);
x0bZeroMean = x0-repmat(permute(mux0,[4 1 2 3]),size(x0,1),size(x0,2));
x1bZeroMean = x0-repmat(permute(mux1,[4 1 2 3]),size(x1,1),size(x1,2));

sz = size(x0);
sz(length(sz)+1:4) = 1;
NCC = zeros([length(grossLags),sz(3:end)]);
x0bZeroMeanReshape = reshape(x0bZeroMean,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]);
for i = 1:length(grossLags);
    idx = max(1,min(K,(1:K)+grossLags(i)));
    x1bZeroMeanReshape = reshape(x1bZeroMean(idx,:,:),[size(x1,1)*size(x1,2),size(x1,3),size(x1,4)]);
    NCC(i,:,:) = nansum(conj(x1bZeroMeanReshape).*x0bZeroMeanReshape)./((sz(1)*sz(2)).*(sigx0.*sigx1));
end
[CCpk CCpkidx] = max(NCC,[],1);
SampleShift = grossLags(CCpkidx);

sigx1SampleShifted = sigx1;
x1SampleShifted = x1;

SampleShifts = repmat(permute(SampleShift,[4 1 2 3]),size(x0,1),size(x0,2));
for j = 1:size(x1,3)
    for k = 1:size(x1,4);
        idx = max(1,min(K,(1:K)+SampleShift(1,j,k)));
        x1SampleShifted(:,:,j,k) = x1(idx,:,j,k).*exp(1j*2*pi*f0*SampleShifts(:,:,j,k)*(1/fs));
    end
end
x1 = x1SampleShifted;

x1Reshape = reshape(x1,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]);
x0Reshape = reshape(x0,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]);
x1ZeroMeanReshape = x1Reshape - repmat(nanmean(x1Reshape,1),size(x1Reshape,1),1);
x0ZeroMeanReshape = x0Reshape - repmat(nanmean(x0Reshape,1),size(x0Reshape,1),1);
corrVal = nansum(conj(x1ZeroMeanReshape).*(x0ZeroMeanReshape))./((sz(1)*sz(2)).*sigx1.*sigx0);
corrValPhasCrct = permute(angle(corrVal),[4 1 2 3]);
tau = corrValPhasCrct/(2*pi*f0);

for n = 1:N-1
    x0PhaseShifted = subsampleshift(x0,f0,tau/2);
    x1PhaseShifted = subsampleshift(x1,f0,-tau/2);
    x1Reshape = reshape(x1PhaseShifted,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]);
    x0Reshape = reshape(x0PhaseShifted,[size(x0,1)*size(x0,2),size(x0,3),size(x0,4)]);
    x1ZeroMeanReshape = x1Reshape - repmat(nanmean(x1Reshape,1),size(x1Reshape,1),1);
    x0ZeroMeanReshape = x0Reshape - repmat(nanmean(x0Reshape,1),size(x0Reshape,1),1);
    corrVal = nansum(conj(x1ZeroMeanReshape).*(x0ZeroMeanReshape),1)./((sz(1)*sz(2)).*sigx1.*sigx0);
    corrValPhasCrct = angle(exp(-1j*2*pi*f0*(tau)/fs).*permute(corrVal,[4 1 2 3]));
    tau = tau + corrValPhasCrct/(2*pi*f0); 
end
if nargout > 1
    CC = squeeze(corrVal);
end

Tau = tau+permute(SampleShift*(1/fs),[4 1 2 3]);
end

