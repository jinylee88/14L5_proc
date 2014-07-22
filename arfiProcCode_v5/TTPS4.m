
function [CT NCC DX] = TTPS4(theta,dtheta,axial,apex,t,sweidata,pushbeamNum,delta,sws);
if ~exist('delta','var')
    delta = 1;
end
[X Z dX] = sectorCoordinates(theta,dtheta,axial,apex);
[vel t1] = differentiateDisplacements(sweidata,t,2000);
vel = removeLateralMedian(vel);
[vL vR] = LeftRightFilter(vel);
DX3 = repmat(permute(dX,[1 2 4 3]),[1 1 size(vel,3) 1]);
vLR = vel;
vLR(DX3<0) = 2*vL(DX3<0);
vLR(DX3>0) = 2*vR(DX3>0);
T3 = repmat(reshape(t1(:),1,1,[]),[size(DX3,1),size(DX3,2),1,size(DX3,4)]);
DXDT = abs(DX3)./T3;
vTest = vLR;
vTest = convn(vTest,ones(5,1,1)/5,'same');
vTest = vTest-repmat(nanmean(vTest,3),[1 1 size(vTest,3)]);
vTest = vTest./repmat(nanstd(vTest,1,3),[1 1 size(vTest,3)]);
vLL = vTest;
vRR = vTest;
vLL(DX3>=-1*sws.xrange(1)) = nan;
vLL(DX3<=-1*sws.xrange(2)) = nan;
vRR(DX3<=1*sws.xrange(1)) = nan;
vRR(DX3>=1*sws.xrange(2)) = nan;
vLL = reshape(vLL,size(vLL,1),size(vLL,2),[]);
vRR = reshape(vRR,size(vRR,1),size(vRR,2),[]);
dX = reshape(dX,size(dX,1),[]);
dt = mean(diff(t));
tLen = size(vRR,3);
nx = size(vRR,2);
Delta = [1];

idxRef = [];
idxTgt = [];
H = [];
for deltaIdx = 1:length(Delta);
    delta = Delta(deltaIdx);
    idx0 = 1:nx-delta;
    idx1 = 1+delta:nx;
    idxRef = [idxRef idx0]; 
    idxTgt = [idxTgt idx1]; 
    for i = 1:length(idx0)
    Hrow = zeros(1,nx-1);
    Hrow(idx0(i):(idx1(i)-1)) = 1;
    H = [H;Hrow];
    end
end
vRefL = vLL(:,idxRef,:);
vTgtL = vLL(:,idxTgt,:);
vTgtL = padarray(vTgtL,[0 0 1],NaN,'both');
vRefR = vRR(:,idxRef,:);
vTgtR = vRR(:,idxTgt,:);
vTgtR = padarray(vTgtR,[0 0 1],NaN,'both');

dxRef = dX(:,idxRef);
dxTgt = dX(:,idxTgt);
swDir = sign(0.5*(dxRef+dxTgt));
dDX = dxTgt-dxRef;

%sz = size(vTgt);
dTtest = 0:0.05:1;
for testIdx = 1:length(dTtest);
    shiftT = dTtest(testIdx);
    shiftN = shiftT/dt;
    shiftI = floor(shiftN);
    shiftF = mod(shiftN,1);
    tidx0 = max(0,min(tLen,(1:tLen)+shiftI))+1;
    tidx1 = max(0,min(tLen,(1:tLen)+shiftI+1))+1;
    vTgtShiftR = ...
        (1-shiftF)*vTgtR(:,:,tidx0)...
        +(shiftF)*vTgtR(:,:,tidx1);
    NCCcoeffR(:,:,testIdx) = nanmean(vRefR.*vTgtShiftR,3);
end
for testIdx = 1:length(dTtest);
    shiftT = -1*dTtest(testIdx);
    shiftN = shiftT/dt;
    shiftI = floor(shiftN);
    shiftF = mod(shiftN,1);
    tidx0 = max(0,min(tLen,(1:tLen)+shiftI))+1;
    tidx1 = max(0,min(tLen,(1:tLen)+shiftI+1))+1;
    vTgtShiftL = ...
        (1-shiftF)*vTgtL(:,:,tidx0)...
        +(shiftF)*vTgtL(:,:,tidx1);
    NCCcoeffL(:,:,testIdx) = nanmean(vRefL.*vTgtShiftL,3);
end
[NCCL(:,:) DTIL(:,:)] = nanmax(NCCcoeffL,[],3);
[NCCR(:,:) DTIR(:,:)] = nanmax(NCCcoeffR,[],3);
DTL(:,:) = dTtest(DTIL(:,:));
DTR(:,:) = dTtest(DTIR(:,:));

%[NCC(:,:,pushIdx) DT(:,:,pushIdx)] = subsamplepeak(dTtest,NCCcoeff,3);
CTL(:,:) = abs(dDX)./DTL(:,:);
CTR(:,:) = abs(dDX)./DTR(:,:);

DX(:,:) = 0.5*(dxRef+dxTgt);;

    

