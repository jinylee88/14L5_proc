function [u, fdemest]=loupas(iref,qref,idisp,qdisp,RF_t,M,fdem,fc,kasai_scale)

% This function was copied from:
% /nefs/ncr2/ProcessArfiData/v3.0/routines/ComputeUpsampledIQdata.m
%
% sjr6 3/12/12

num1=zeros(RF_t,1);
den1=zeros(RF_t,1);
num2=zeros(RF_t,1);
den2=zeros(RF_t,1);
tmp =zeros(RF_t,1);

for k=1:RF_t-1, tmp(k)=qref(k)*iref(k+1)-iref(k)*qref(k+1); end
num1=mom1v(tmp,RF_t-1,M-1,num1);
for k=1:RF_t-1, tmp(k)=qdisp(k)*idisp(k+1)-idisp(k)*qdisp(k+1); end
num2=mom1v(tmp,RF_t-1,M-1,num2);
for k=1:RF_t-M, num2(k)=num2(k)+num1(k); end

for k=1:RF_t-1, tmp(k)=iref(k)*iref(k+1)+qref(k)*qref(k+1); end
den1=mom1v(tmp,RF_t-1,M-1,den1);
for k=1:RF_t-1, tmp(k)=idisp(k)*idisp(k+1)+qdisp(k)*qdisp(k+1); end
den2=mom1v(tmp,RF_t-1,M-1,den2);
for k=1:RF_t-M, den2(k)=den2(k)+den1(k); end

for k=1:RF_t, tmp(k)=qref(k)*idisp(k)-iref(k)*qdisp(k); end
num1=mom1v(tmp,RF_t,M,num1);
for k=1:RF_t, tmp(k)=iref(k)*idisp(k)+qref(k)*qdisp(k); end
den1=mom1v(tmp,RF_t,M,den1);

fdem = reshape(fdem, size(num1));
fc = reshape(fc, size(num1));
u=atan2(num1,den1)./(1+(1./fdem).*atan2(num2,den2)/(2*pi))*kasai_scale./fc./2;

if nargout==2,fdemest = atan2(num2,den2)/(2*pi);end