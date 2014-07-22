function [Tau CC]= PesaventoRF(x1b,x2b,fs,f0,Klen,N)
if ~exist('Klen','var')
    Klen = 1;
end
if ~exist('N','var')
    N = 3;
end
x1b = repmat(x1b,size(x2b)./size(x1b));
tau = zeros([1 size(x1b)]);
Tau = zeros(size(x1b));
K = size(x1b,1);
if nargout > 1
    CC = zeros(size(x1b));
end
x1b1 = reshape(x1b,[1 size(x1b)]);
x2b1 = reshape(x2b,[1 size(x1b)]);
lags = -floor(Klen/2):floor(Klen/2);
Klen = length(lags);
for i = 1:length(lags);
    x1b1(i,:,:,:) = x1b(max(1,min(K,(1:K)+lags(i))),:,:);
    x2b1(i,:,:,:) = x2b(max(1,min(K,(1:K)+lags(i))),:,:);
end
x1bt = x1b1;
x2bt = x2b1;
if 1 || Klen == 1
    corrVal = sum(conj(x2bt).*x1bt,1);
else
    corrVal = sum(conj(x2bt-repmat(mean(x2bt,1),[Klen 1 1 1])).*(x1bt-repmat(mean(x1bt,1),[Klen 1 1 1])))./((Klen).*(std(conj(x2bt),1,1).*std(x1bt,1,1)));
end
corrValPhasCrct = angle(exp(0*-1j*2*pi*f0*(tau)/fs).*corrVal);
tau = tau + corrValPhasCrct/(2*pi*f0);
cc(:,1) = corrVal;
T(:,1) = tau;
for n = 1:N-1
    x1bt = subsampleshift(x1b1,f0,tau/2);
    x2bt = subsampleshift(x2b1,f0,-tau/2);
    corrVal0 = corrVal;
    if 1 || Klen == 1
        corrVal = sum(conj(x2bt).*x1bt,1);
    else
        corrVal = sum(conj(x2bt-repmat(mean(x2bt,1),[Klen 1 1 1])).*(x1bt-repmat(mean(x1bt,1),[Klen 1 1 1])))./((Klen).*(std(conj(x2bt),1,1).*std(x1bt,1,1)));
    end
    corrValPhasCrct = angle(exp(0*1j*2*pi*f0*(tau)/fs).*corrVal);
    tau = tau + corrValPhasCrct/(2*pi*f0); 
    cc(:,1+n) = corrVal;
    T(:,1+n) = tau;
end
if nargout > 1
    CC = squeeze(corrVal);
end

Tau = squeeze(tau);
end

