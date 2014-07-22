
function [CT NCC DX] = TTPS3(theta,dtheta,axial,apex,t,sweidata,pushbeamNum,delta,sws);
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
%vTest(DXDT<sws.dxdtrange(1)) = nan;
%vTest(DXDT>sws.dxdtrange(2)) = nan;
vTest(abs(DX3)<sws.xrange(1)) = nan;
vTest(abs(DX3)>sws.xrange(2)) = nan;
%vTest(T3<sws.trange(1)) = nan;
%vTest(T3>sws.trange(2)) = nan;

dt = mean(diff(t));
cTtest = 0.5:0.1:8;
tIdxRef = find(t1>-inf & t1<inf);
k0 = tIdxRef(1);
tLen = length(tIdxRef);
for pushIdx = 1:size(vTest,4)
idxR0 = ceil(pushbeamNum(pushIdx)):1:size(dX,2)-delta;
idxR1 = ceil(pushbeamNum(pushIdx))+delta:1:size(dX,2);
idxL0 = floor(pushbeamNum(pushIdx)):-1:1+delta;
idxL1 = floor(pushbeamNum(pushIdx))-delta:-1:1;


firstRight = find((1:size(dX,2))>=pushbeamNum(pushIdx),1,'first');
idx1 = 1:size(dX,2);
idx0 = [min(floor(pushbeamNum(pushIdx)),(1:firstRight-1)+delta) max(ceil(pushbeamNum(pushIdx)),(firstRight:size(dX,2))-delta)];
%idx0 = [fliplr(idxL0) idxR0];
%idx1 = [fliplr(idxL1) idxR1];
%[I J K0] = ndgrid(1:size(vTest,1),1:length(idx0),tIdxRef);
vRef = vTest(:,idx0,:,pushIdx);
vTgt = vTest(:,idx1,:,pushIdx);
vTgt = padarray(vTgt,[0 0 1],NaN,'both');
dxRef = dX(:,idx0,pushIdx);
dxTgt = dX(:,idx1,pushIdx);
swDir = sign(0.5*(dxRef+dxTgt));
dDX = dxTgt-dxRef;
cTtest = 0.5:0.1:10;
%sz = size(vTgt);
dTtest = -0.5:0.01:0.5;
for testIdx = 1:length(dTtest);
    %cT = cTtest(testIdx);
    %shiftT = abs(dDX)./cT;
    shiftT = dTtest(testIdx);
    shiftN = shiftT/dt;
    shiftI = floor(shiftN);
    shiftF = mod(shiftN,1);
    tidx0 = max(0,min(tLen,(1:tLen)+shiftI))+1;
    tidx1 = max(0,min(tLen,(1:tLen)+shiftI+1))+1;
    vTgtShift = ...
        (1-shiftF)*vTgt(:,:,tidx0)...
        +(shiftF)*vTgt(:,:,tidx1);
    NCCcoeff(:,:,testIdx) = nanmean(vRef.*vTgtShift,3);
%  subplot(211)
%  imagesc(squeeze(vRef(80,:,:)),[-10 10]);
%  subplot(212)
%  imagesc(squeeze(vTgtShift(80,:,:)),[-10 10]);
% %set(im1,'cdata',squeeze(vRef1(30,:,:))');
% %set(im2,'cdata',squeeze(vTgtShift(30,:,:))');
%  pause(0.01);
end
%[NCC(:,:,pushIdx) DTI(:,:,pushIdx)] = nanmax(NCCcoeff,[],3);
%DT(:,:,pushIdx) = dTtest(DTI(:,:,pushIdx));
[NCC(:,:,pushIdx) DT(:,:,pushIdx)] = subsamplepeak(dTtest,convn(NCCcoeff,ones(1,3,1),'same'),3);
CT(:,:,pushIdx) = abs(dDX)./DT(:,:,pushIdx);
DX(:,:,pushIdx) = 0.5*(dxRef+dxTgt);;
end
K = ones(sws.smoothingkernel(1),sws.smoothingkernel(2),1);
CT1 = CT;
CT1(isnan(CT)) = 0;
CT1 = convn(CT1.*NCC,K,'same')./convn(~isnan(CT),K,'same');
CT = CT1;
fprintf('done')

