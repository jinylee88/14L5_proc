function [cc_coef] = computeCC(cData,kernelLength,refidx)

% sjr6 adapted from pjh7 on 5/16/12

if nargin<3
    refidx = 1;
end
k = ones(kernelLength,1);
nt = size(cData,3);
ref = repmat(cData(:,:,refidx),[1 1 nt]);
cc_coef = convn(conj(ref).*cData,k,'same')./sqrt(repmat(convn(abs(ref(:,:,1)).^2,k,'same'),[1 1 nt]).*convn(abs(cData).^2,k,'same'));
cc_coef = cc_coef(1:size(cc_coef,1)-(kernelLength+1),:,:);

