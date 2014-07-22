function [Tau CC]= PesaventoParallel(x1b,x2b,fs,f0,Klen,N)
if ~exist('Klen','var')
Klen = 1;
end
if ~exist('N','var')
N = 3;
end
tau = zeros(1,size(x1b,2),size(x1b,3));
Tau = zeros(size(x1b));
K = size(x1b,1);
if nargout > 1
CC = zeros(size(x1b));
end
for k = 1:K;
    kidx = max(1,k-floor(Klen/2)):min(K,k+floor(Klen/2));
    for n = 1:N
    x1bt = subsampleshift(x1b(kidx,:,:),f0,tau(ones(length(kidx),1),:,:)/2);
    x2bt = subsampleshift(x2b(kidx,:,:),f0,-tau(ones(length(kidx),1),:,:)/2);
    if length(kidx) == 1
    corrVal = sum(conj(x2bt).*x1bt,1);
    else
    corrVal = sum(conj(x2bt-repmat(mean(x2bt,1),[length(kidx) 1 1])).*(x1bt-repmat(mean(x1bt,1),[length(kidx) 1 1])))./((length(kidx)).*(std(conj(x2bt),1,1).*std(x1bt,1,1)));
    end
    %corrVal = ncorrAccel_V3(conj(x2bt),x1bt);
    corrValPhasCrct = angle(exp(-1j*2*pi*f0*(tau)/fs).*corrVal);
    tau = tau + corrValPhasCrct/(2*pi*f0);
    end
     if nargout > 1    
    x1bt = subsampleshift(x1b(kidx,:,:),f0,tau(ones(length(kidx),1),:,:)/2);
    x2bt = subsampleshift(x2b(kidx,:,:),f0,-tau(ones(length(kidx),1),:,:)/2);
    %CC(k,:,:) = ncorrAccel_V3(conj(x2bt),x1bt);
     if length(kidx) == 1
    corrVal = sum(conj(x2bt).*x1bt,1);
    else
    CC(k,:,:) = sum(conj(x2bt-repmat(mean(x2bt,1),[length(kidx) 1 1])).*(x1bt-repmat(mean(x1bt,1),[length(kidx) 1 1])))./((length(kidx)).*(std(conj(x2bt),1,1).*std(x1bt,1,1)));
     end
     end
    Tau(k,:,:) = tau;
end
end

function x = subsampleshift(x0,f0,tau)
if size(x0,1) == size(tau,1);
x = (ifft(fft(x0).*exp(-1j*2*pi*f0*tau)));
else
x = (ifft(fft(x0)*exp(-1j*2*pi*f0*tau)));
end
end

