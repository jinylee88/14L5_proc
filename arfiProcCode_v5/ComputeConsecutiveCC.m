function [cc_coef] = ComputeConsecutiveCC(data,k_length)
k = ones(k_length,1);
nt = size(data,3);
ref = circshift(data,[0 0 1]);
cc_coef = convn(conj(ref).*data,k,'same')./sqrt(repmat(convn(abs(ref(:,:,1)).^2,k,'same'),[1 1 nt]).*convn(abs(data).^2,k,'same'));
cc_coef(:,:,1) = 1;

