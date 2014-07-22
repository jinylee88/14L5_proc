
function CT = sws_DTTPS(TTP,DX,n)
if ~exist('n','var')
    n = 1;
end
B = [-1*ones(1,ceil(n/2)) zeros(1,mod(n+1,2)) ones(1,ceil(n/2))];
CT = convn(DX,B,'same')./convn(sign(DX).*TTP,B,'same');
end
