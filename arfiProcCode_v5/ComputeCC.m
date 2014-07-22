function [cc_coef] = ComputeCC(data,k_length,refidx)
if ~exist('num_pretracks','var')
    refidx = 1;
end
k = ones(k_length,1);
nt = size(data,3);
ref = repmat(data(:,:,refidx),[1 1 nt]);
for i = 1:100:nt
    idx = i:min(i+99,nt);
    idx = idx(idx<=nt);
cc_coef(:,:,idx) = convn(conj(ref(:,:,idx)).*data(:,:,idx),k,'same')./sqrt(repmat(convn(abs(ref(:,:,1)).^2,k,'same'),[1 1 length(idx)]).*convn(abs(data(:,:,idx)).^2,k,'same'));
end
