function out=upsampleSplineSevalMatlab(n,x,y,xx)

%
% this routine is a combination of the spline and seval functions
%   in GFP's Loupas C code, translated to Matlab.
%
%    WARNING - WARNING - WARNING - WARNING
%
% the seval portion of this routine only works for the variable
%   'x' with integer spacing, i.e., x = 0, 1, 2, 3, ...
%
% do not use this function for general splines
%

% This function was copied from:
% /nefs/ncr2/ProcessArfiData/v3.0/routines/UpsampleSplineSevalMatlab.m
%
% sjr6 3/12/12

if n<=3                         % copy of GFP's 'spline.f' function
    error('n is too small')     %   translated to Matlab
end                             %   see /data/jjd/loupas/spline.f

b=zeros(1,n);
c=zeros(1,n);
d=zeros(1,n);

nm1=n-1;

d(1)=x(2)-x(1);
c(2)=(y(2)-y(1))/d(1);
for i=2:nm1                 % do 10
     d(i)=x(i+1)-x(i);
     b(i)=2.*(d(i-1)+d(i));
     c(i+1)=(y(i+1)-y(i))/d(i);
     c(i)=c(i+1)-c(i);
end

 b(1) = -d(1);
  b(n) = -d(n-1);
  c(1) = 0.;
  c(n) = 0.;
%      if ( n .eq. 3 ) go to 15;
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1));
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3));
  c(1) = c(1)*d(1)^2/(x(4)-x(1));
  c(n) = -c(n)*d(n-1)^2/(x(n)-x(n-3));

%   do 20 i = 2, n
for i = 2:n
     t = d(i-1)/b(i-1);
     b(i) = b(i) - t*d(i-1);
     c(i) = c(i) - t*c(i-1);
end

c(n) = c(n)/b(n);
%      do 30 ib = 1, nm1
for ib = 1:nm1
     i = n-ib;
     c(i) = (c(i) - d(i)*c(i+1))/b(i);
end

b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n));
%      do 40 i = 1, nm1
for i = 1:nm1
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i));
     d(i) = (c(i+1) - c(i))/d(i);
     c(i) = 3.*c(i);
end
  c(n) = 3.*c(n);
  d(n) = d(n-1);

    %
    % end of spline.f function, now to seval.f function
    %
    %   this has been changed compared to /data/jjd/loupas/seval.f
    %   to be optimized for integer spacing, i.e., x=0,1,2,3,...
    %   
    % WARNING - WARNING - WARNING
    %
    %   do not use this routine for general spline evaluation
    %

out=zeros(length(xx),1);

for ixx=1:length(xx)
    xval=floor(xx(ixx));
    i=xval+1;               % get index for y,b,c,d based on x value
    dx=xx(ixx)-xval;        
    out(ixx)=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)));    
end
