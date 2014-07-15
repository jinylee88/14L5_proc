function [u]=loupasParallel(iref,qref,idisp,qdisp,RF_t,M,fdem,fc,kasai_scale)

% This function was copied from:
% /nefs/ncr2/ProcessArfiData/v3.0/routines/ComputeUpsampledIQdata.m
%
% sjr6 3/12/12
%
% This was adapted to work with a parfor loop
% This was also modified to use vector math rather than for loops

u = zeros(size(idisp));
for i = 1:size(idisp,2)
    idisptmp = idisp(:,i);
    qdisptmp = qdisp(:,i);
    
    filtKern1 = ones(1,M-1)/(M-1);
    filtKern2 = ones(1,M)/(M);
    
    tmp = qref(1:RF_t-1).*iref(2:RF_t) - iref(1:RF_t-1).*qref(2:RF_t);
    num1 = filter(filtKern1, 1, tmp);
    
    tmp = qdisptmp(1:RF_t-1).*idisptmp(2:RF_t) - idisptmp(1:RF_t-1).*qdisptmp(2:RF_t);
    num2 = filter(filtKern1, 1, tmp);
    num2 = num2+num1;
    
    tmp = iref(1:RF_t-1).*iref(2:RF_t) + qref(1:RF_t-1).*qref(2:RF_t);
    den1 = filter(filtKern1, 1, tmp);
    
    tmp = idisptmp(1:RF_t-1).*idisptmp(2:RF_t) + qdisptmp(1:RF_t-1).*qdisptmp(2:RF_t);
    den2 = filter(filtKern1, 1, tmp);
    den2 = den2+den1;
    
    tmp = qref.*idisptmp-iref.*qdisptmp;
    num1 = filter(filtKern2, 1, tmp);
    
    tmp = iref.*idisptmp+qref.*qdisptmp;
    den1 = filter(filtKern2, 1, tmp);
    
    vec = 1:RF_t-M;
    u(vec,i) = fdem(vec).*atan2(num1(vec+M-1),den1(vec+M-1))...
                        ./ (fdem(vec) + atan2(num2(vec+M-2),den2(vec+M-2))/(2*pi))...
                        .* kasai_scale./fc(vec)./2;

end
% for i = 1:size(idisp,2)
%     idisptmp = idisp(:,i);
%     qdisptmp = qdisp(:,i);
% 
%     num1=zeros(RF_t,1);
%     den1=zeros(RF_t,1);
%     num2=zeros(RF_t,1);
%     den2=zeros(RF_t,1);
%     tmp =zeros(RF_t,1);
% 
%     for k=1:RF_t-1, tmp(k)=qref(k)*iref(k+1)-iref(k)*qref(k+1); end
%     num1=mom1v(tmp,RF_t-1,M-1,num1);
%     for k=1:RF_t-1, tmp(k)=qdisptmp(k)*idisptmp(k+1)-idisptmp(k)*qdisptmp(k+1); end
%     num2=mom1v(tmp,RF_t-1,M-1,num2);
%     for k=1:RF_t-M, num2(k)=num2(k)+num1(k); end
% 
%     for k=1:RF_t-1, tmp(k)=iref(k)*iref(k+1)+qref(k)*qref(k+1); end
%     den1=mom1v(tmp,RF_t-1,M-1,den1);
%     for k=1:RF_t-1, tmp(k)=idisptmp(k)*idisptmp(k+1)+qdisptmp(k)*qdisptmp(k+1); end
%     den2=mom1v(tmp,RF_t-1,M-1,den2);
%     for k=1:RF_t-M, den2(k)=den2(k)+den1(k); end
% 
%     for k=1:RF_t, tmp(k)=qref(k)*idisptmp(k)-iref(k)*qdisptmp(k); end
%     num1=mom1v(tmp,RF_t,M,num1);
%     for k=1:RF_t, tmp(k)=iref(k)*idisptmp(k)+qref(k)*qdisptmp(k); end
%     den1=mom1v(tmp,RF_t,M,den1);
% 
%     fdem = reshape(fdem, size(num1));
%     fc = reshape(fc, size(num1));
%     u(:,i)=atan2(num1,den1)./(1+(1./fdem).*atan2(num2,den2)/(2*pi))*kasai_scale./fc./2;
% 
% end