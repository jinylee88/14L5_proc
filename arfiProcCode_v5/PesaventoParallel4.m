function [Tau CC]= PesaventoParallel4(IQ0,IQ1,fs,f0,Klen,SrchWin,N)
if ~exist('Klen','var')
    Klen = 1;
end
if ~exist('SrchWin','var')
    SrchWin = 100e-6/770*fs;
end
if ~exist('N','var')
    N = 3;
end
grossLags = [-ceil(SrchWin):ceil(SrchWin)];
IQ0 = repmat(IQ0,size(IQ1)./size(IQ0));
tau = zeros([1 size(IQ0)]);
Tau = zeros(size(IQ0));
K = size(IQ0,1);
if nargout > 1
    CC = zeros(size(IQ0));
end
kernelLags = -floor(Klen/2):floor(Klen/2);
Klen = length(kernelLags);
tic
fprintf('Axial Sample %03.0f/%03.0f',0,size(IQ0,1))
parfor axidx = 1:size(IQ0,1);
    fprintf('%s%03.0f/%03.0f',8*ones(1,7),axidx,size(IQ0,1));
    idx = max(1,axidx+kernelLags(1)):min(K,axidx+kernelLags(end));
    Klen = length(idx);
    x0 = IQ0(idx,:,:);
    mux0 = nanmean(x0,1);
    sigx0 = nanstd(x0,1,1);
    x0bZeroMean = x0-mux0(ones(Klen,1),:,:);
    x0bZeroMean = x0bZeroMean./sigx0(ones(Klen,1),:,:);
    
%     if length(grossLags)==1 && grossLags ==0;
        x1 = IQ1(idx,:,:);
        mux1 = nanmean(x1,1);
        sigx1 = nanstd(x1,1,1);
        x1bZeroMean = x1-mux1(ones(Klen,1),:,:,:);
        x1bZeroMean = x1bZeroMean./sigx1(ones(Klen,1),:,:);
        SampleShift = 0;
%     else
%         
%         for grossLagIdx = 1:length(grossLags)
%             x1 = IQ1(max(1,min(K,idx+grossLags(grossLagIdx))),:,:);
%             mux1 = nanmean(x1,1);
%             sigx1 = nanstd(conj(x1),1,1);
%             x1bZeroMean = x1-mux1(ones(Klen,1),:,:,:);
%             NCCcoeff(grossLagIdx,:,:) = nansum(conj(x1bZeroMean.*x0bZeroMean))./((Klen).*(sigx0.*sigx1));
%         end
%         [NCC grossLagPeakIdx] = max(abs(NCCcoeff),[],1);
%         SampleShift = grossLags(grossLagPeakIdx);
%         
%         x1SampleShifted = x1;
%         for j = 1:size(x1,2)
%             for k = 1:size(x1,3);
%                 idx1 = max(1,min(K,idx+SampleShift(1,j,k)));
%                 x1SampleShifted(:,j,k) = IQ1(idx1,j,k).*exp(1j*2*pi*f0*SampleShift(1,j,k)*(1/fs));
%             end
%         end
%         x1 = x1SampleShifted;
%         sigx1 = nanstd(x1SampleShifted,1,1);
%         x1bZeroMean = x1-mux1(ones(Klen,1),:,:,:);
%         x1bZeroMean = x1bZeroMean./sigx1(ones(Klen,1),:,:);
%     end
    
    corrVal = nanmean(conj(x1bZeroMean).*(x0bZeroMean),1);
    corrValPhasCrct = angle(corrVal);
    tau = corrValPhasCrct/(2*pi*f0);

    for n = 1:N-1
        x0PhaseShifted = subsampleshift(x0bZeroMean,f0,tau/2);
        x1PhaseShifted = subsampleshift(x1bZeroMean,f0,-tau/2);
        %corrVal = nansum(conj(x1PhaseShifted-repmat(nanmean(x1PhaseShifted,1),[Klen 1 1 1])).*(x0PhaseShifted-repmat(nanmean(x0PhaseShifted,1),[Klen 1 1 1])))./((Klen).*sigx1.*sigx0);
        corrVal = nanmean(conj(x1PhaseShifted).*(x0PhaseShifted),1);
        corrValPhasCrct = angle(exp(-1j*2*pi*f0*(tau)/fs).*corrVal);
        tau = tau + corrValPhasCrct/(2*pi*f0);
    end
    CC(axidx,:,:) = corrVal;
    Tau(axidx,:,:) = tau+(SampleShift*(1/fs));
end
fprintf('...done (%0.1fs)\n',toc);
