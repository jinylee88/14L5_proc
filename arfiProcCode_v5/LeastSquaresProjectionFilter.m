function B = LeastSquaresProjectionFilter(K,gamma)
if ~exist('gamma','var')
    gamma = 0;
end
if length(K) == 2;
    k0 = K(1);
    k = K(2);
else
    k0 = 1;
    k = K;
end
if k<k0
    warning('input:kltk','k must be in ascending order')
    k = k0;
end
A = zeros(1,k);
ii = 1;
for i = 1:k
    for j = i+(k0-1):k
        A(ii,i:j) = 1;
        ii = ii+1;
    end
end
L  = 0.5*(eye(k-1,k) - circshift(eye(k-1,k),[0 1]));
coeff = (inv(A.'*A+gamma.*(L'*L))*(A.'));
%coeff = (inv(A.'*A+gamma.*eye(k))*(A.'));

B = zeros(k+1,k+1);
ii = 1;
for i = 1:k
    for j = i+(k0-1):k
        B(i,j+1) = coeff((k+1)/2,ii);
        B(j+1,i) = -1*coeff((k+1)/2,ii);
        ii = ii+1;
    end
end
%B = B./sum(B(:));
