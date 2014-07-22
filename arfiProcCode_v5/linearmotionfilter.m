function [arfidata motion A] = linearmotionfilter(arfidata0,t,tidx,order,cc,ccthresh)
arfidata0 = 1e6*arfidata0;
arfidata1 = arfidata0;
if exist('ccthresh','var');
mask = cc<ccthresh;
arfidata1(arfidata1==0) = 1e-15;
arfidata1(mask) = 0;
while any(mask)
arfidata2 = convn(arfidata1,ones(3,1,3),'same');
counts = convn(double(~mask),ones(3,1,3),'same');
arfidata1(mask) = arfidata2(mask(:))./counts(mask(:));
arfidata1((abs(arfidata1)==inf)|isnan(arfidata1)) = 0;
mask = (arfidata1 == 0);
end

end
t = t./mode(diff(t));
%arfidata1 = 1e6 * arfidata1;
if order == 1 %don't bother inverting the matrix, we can skip to the least squares solution
Y = arfidata1(:,:,tidx);
X = repmat(reshape(t(tidx),1,1,[]),[size(Y,1),size(Y,2),1]);
muY = repmat(mean(Y,3),[1,1,length(tidx)]);
muX = repmat(mean(X,3),[1,1,length(tidx)]);
Beta = repmat(sum((X-muX).*(Y-muY),3)./sum((X-muX).^2,3),[1 1 length(t)]);
Alpha = repmat(muY(:,:,1)-Beta(:,:,1).*muX(:,:,1),[ 1 1 length(t)]);
T = repmat(permute(t,[1 3 2]),[size(Y,1) size(Y,2) 1]);
motion = Alpha + Beta.*T;
arfidata = 1e-6*(arfidata0-motion);
motion = 1e-6*motion;
A(:,:,1) = Alpha(:,:,1);
A(:,:,2) = Beta(:,:,1);

else
    
Y = arfidata1(:,:,tidx); 
Y = permute(Y,[3 1 2]);   
Y = reshape(Y,length(tidx),[]);
    X = [];
    T = [];
for i = 0:order
    X = [X t(tidx).^i'];
    T = [T t.^i'];
end
XtXiXt = inv(X'*X)*X';

A = XtXiXt*Y;
motion = T*A;
motion = reshape(motion,length(t),size(arfidata0,1),size(arfidata0,2));
motion = permute(motion,[2 3 1]);

% tidx1 = find(t<-1.5);
% Y = arfidata1(:,:,tidx1);
% X = repmat(reshape(t(tidx1),1,1,[]),[size(Y,1),size(Y,2),1]);
% muY = repmat(mean(Y,3),[1,1,length(tidx1)]);
% muX = repmat(mean(X,3),[1,1,length(tidx1)]);
% Beta = repmat(sum((X-muX).*(Y-muY),3)./sum((X-muX).^2,3),[1 1 length(t)]);
% Alpha = repmat(muY(:,:,1)-Beta(:,:,1).*muX(:,:,1),[ 1 1 length(t)]);
% T = repmat(permute(t,[1 3 2]),[size(Y,1) size(Y,2) 1]);
% linmotion = Alpha + Beta.*T;
% 
% w = max(0,min(1,t./(-2*min(t))));
% W = repmat(permute(w,[1 3 2]),size(arfidata0,1),size(arfidata0,2) );
% cmotion = W.*motion + (1-W).*linmotion;


arfidata = 1e-6*(arfidata0-motion);
motion = 1e-6*motion;
end