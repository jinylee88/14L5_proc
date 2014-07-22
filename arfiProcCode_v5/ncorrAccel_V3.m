function nCorr = ncorrAccel_V3(ker,srch);
%nCorr = ncorrAccel_V3(ker,srch) Normalized Cross Correlation
%
%Finds Normalized Cross-Correlation Coefficient in the first dimension.
%is size(srch,1) > size(ker,1), the coefficients will be calculated for all
%possible lags in the first dimension.

n = size(ker,1);
kerMean = mean(ker,1);
srchMean = convn(srch,ones(n,1,1),'valid')/n;
kerEnergy = std(ker,1,1);
srchEnergy = sqrt(convn(srch.^2,ones(n,1,1),'valid')/n-srchMean.^2);
N = size(srchMean,1);
nCorr = zeros(size(srchMean));
for i = 1:N
    ker0 = ker-repmat(kerMean,[n,1,1]);
    srch0 = srch(i+(0:n-1),:,:)-repmat(srchMean(i,:,:),[n,1,1]);
    nCorr(i,:,:) = (1/n) * sum(ker0.*srch0,1)./(kerEnergy.*srchEnergy(i,:,:));
end


