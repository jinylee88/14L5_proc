function [iUp,qUp]=computeUpsampledIQdata(iz,qz,upsampleFactor)

%
% this function calculates the upsampled iz,qz data using the
% factor upsampleFactor
%

% This function was copied from:
% /nefs/ncr2/ProcessArfiData/v3.0/routines/ComputeUpsampledIQdata.m
%
% sjr6 3/12/12

m=numel(iz);       % get number of elements, OK or 1D or 3D
n=m*upsampleFactor;

x=(0:m-1);                      % xvals of iz,qz data
xx=(0:n-1)/upsampleFactor;      % xvals of iUp,qUp data

iz=reshape(iz,m,1);             % reshape to vector
qz=reshape(qz,m,1);

iUp=upsampleSplineSevalMatlab(m,x,iz,xx);
qUp=upsampleSplineSevalMatlab(m,x,qz,xx);
