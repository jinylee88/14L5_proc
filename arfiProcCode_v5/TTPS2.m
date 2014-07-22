
function [TTP DX vTest] = TTPS2(theta,dtheta,axial,apex,t,sweidata,sws);
[X Z DX] = sectorCoordinates(theta,dtheta,axial,apex);
[vel t1] = differentiateDisplacements(sweidata,t,2000);
vel = removeLateralMedian(vel);
[vL vR] = LeftRightFilter(vel);
DX3 = repmat(permute(DX,[1 2 4 3]),[1 1 size(vel,3) 1]);
vLR = vel;
vLR(DX3<0) = 2*vL(DX3<0);
vLR(DX3>0) = 2*vR(DX3>0);
T3 = repmat(reshape(t1(:),1,1,[]),[size(DX3,1),size(DX3,2),1,size(DX3,4)]);
DXDT = abs(DX3)./T3;
vTest = vLR;
vTest = convn(vTest,ones(5,1,1)/5,'same');
dt = mean(diff(t));
cTtest = 0.5:0.1:8;
tIdxRef = find(t1>-inf & t1<inf);
k0 = tIdxRef(1);
tLen = length(tIdxRef);
[I J K0] = ndgrid(1:size(vTest,1),1:size(vTest,2)-1,tIdxRef);
for pushIdx = 1:size(vTest,4)
vRef = vTest(:,1:end-1,tIdxRef,pushIdx);
vTgt = vTest(:,2:end,:,pushIdx);
dxRef = DX(:,1:end-1,pushIdx);
dxTgt = DX(:,2:end,pushIdx);
swDir = sign(0.5*(dxRef+dxTgt));
dDX = diff(DX(:,:,pushIdx),1,2);
maxLag = max(dxRef(:)-dxTgt(:));

sz = size(vTgt);
dTtest = linspace(0,0.25,length(cTtest));
for testIdx = 1:length(cTtest);
    cT = cTtest(testIdx);
    shiftT = swDir.*dDX./cT;
    shiftT = swDir*dTtest(testIdx);
    shiftN = shiftT/dt;
    shiftI = floor(shiftN);
    shiftF = mod(shiftN,1);
    K1 = K0 + repmat(shiftI,[1 1 tLen]);
    K2 = K1+1;
    Kvalid = K1<=(sz(3)) & K1>=1 & K2<=(sz(3)) & K2>=1;
    vTgtShift1 = vTgt(sub2ind(sz,I,J,max(1,min(sz(3),K1))));
    vTgtShift2 = vTgt(sub2ind(sz,I,J,max(1,min(sz(3),K2))));
    vTgtShift = vTgtShift1.*repmat((1-shiftF),[1 1 tLen]) + vTgtShift2.*repmat((shiftF),[1 1 tLen]);
    vTgtShift(~Kvalid) = nan;
    vRef1 = vRef;
    vRef1(~Kvalid) = nan;
    muRef = nanmean(vRef1,3);
    muTgt = nanmean(vTgtShift,3);
    sigRef = nanstd(vRef1,1,3);
    sigTgt = nanstd(vTgtShift,1,3);
    NCCcoeff(:,:,testIdx) = (1./sum(Kvalid,3)) .* nansum((vRef1-muRef(:,:,ones(1,tLen))).*(vTgtShift-muTgt(:,:,ones(1,tLen))),3)./(sigRef.*sigTgt);
    CCcoeff(:,:,testIdx) = (1./sum(Kvalid,3)) .* nansum(vRef.*vTgtShift,3);
    
    %CC(:,:,testIdx) = permute(ncorrAccel_V3(permute(vRef,[3 1 2]),permute(vTgtShift,[3 1 2])),[2 3 1]);
%     for i = 1:31
%         subplot(ax(i));
%         plot(1:tLen,squeeze(vRef(30,i,:)),'b-',1:tLen,squeeze(vTgtShift(30,i,:)),'r-');
%         axis tight
%         set(ax(i),'XTick',[],'YTick',[])
%     end
%subplot(211)
%imagesc(squeeze(vRef1(30,:,:)),[-10 10]);
%subplot(212)
%imagesc(squeeze(vTgtShift(30,:,:)),[-10 10]);
set(im1,'cdata',squeeze(vRef1(30,:,:))');
set(im2,'cdata',squeeze(vTgtShift(30,:,:))');
pause(0.1);
end
[CC(:,:,pushIdx) CT(:,:,pushIdx)] = subsamplepeak(cTtest,CCcoeff,3);
[NCC(:,:,pushIdx) CT(:,:,pushIdx)] = subsamplepeak(cTtest,NCCcoeff,3);
[NCC(:,:,pushIdx) DT(:,:,pushIdx)] = subsamplepeak(dTtest,NCCcoeff,3);

end
    


vTest(DXDT<sws.dxdtrange(1)) = nan;
vTest(DXDT>sws.dxdtrange(2)) = nan;
vTest(abs(DX3)<sws.xrange(1)) = nan;
vTest(abs(DX3)>sws.xrange(2)) = nan;
vTest(T3<sws.trange(1)) = nan;
vTest(T3>sws.trange(2)) = nan;
[pk TTP0] = subsamplepeak(t1,vTest,3);
TTP0(pk<(0.75*nanmedian(pk(:)))) = nan;
TTP0(pk>(10*nanmedian(pk(:)))) = nan;
TTP0 = squeeze(TTP0);
for pushIdx = 1:size(TTP0,3);
    TTP(:,:,pushIdx) = mediannan(TTP0(:,:,pushIdx),[sws.smoothingkernel]);
end
end
end

function B = subMatrix(A,k0,kLen)
   B = A(sub2ind(size(A),repmat((1:size(A,1))',[1 size(A,2) size(K,3)]),repmat(1:size(A,2),[size(A,1),1,size(K,3)]),K));
end



