
function [TTP DX vLR] = TTPS(theta,dtheta,axial,apex,t,sweidata,sws);
[X Z DX] = sectorCoordinates(theta,dtheta,axial,apex);
[vel t1] = differentiateDisplacements(sweidata,t,2000);
vel = removeLateralMedian(vel);
 [vL vR] = LeftRightFilter(vel);
 DX3 = repmat(permute(DX,[1 2 4 3]),[1 1 size(vel,3) 1]);
 vLR = vel;
 vLR(DX3<0) = 2*vL(DX3<0);
 vLR(DX3>0) = 2*vR(DX3>0);
 T3 = repmat(reshape(t1(:),1,1,[]),[size(DX3,1),size(DX3,2),1,size(DX3,4)]);
vTest = vLR;
vTest(abs(DX3)<sws.xRange(1)) = nan;
vTest(abs(DX3)>sws.xRange(2)) = nan;
vTest(T3<sws.tRange(1)) = nan;
vTest(T3>sws.tRange(2)) = nan;
[pk TTP0] = subsamplepeak(t1,vTest,3);
TTP0 = squeeze(TTP0);
for i = 1:size(TTP0,3);
    TTP(:,:,i) = mediannan(TTP0(:,:,i),[sws.kernel]);
end
end
