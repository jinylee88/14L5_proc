function [DT CC DX Q] = computeTimeDelay_mex(data,dxdtRange,kerLen)
addpath /getlab/pjh7/SC2000/arfiProcCode_v5/ncorr2
xt = data.lat;
xp = data.lat1;
t = data.t;
z = data.axial;
numtrackloc = length(xt);
numpushloc = length(xp);
CC = single(zeros(length(z),numtrackloc,numtrackloc,numpushloc));
DT = single(zeros(length(z),numtrackloc,numtrackloc,numpushloc));
if nargout>3
    Q = single(zeros(length(z),numtrackloc,numtrackloc,pushloc,3));
    saveQ = 1;
else
    saveQ = 0;
end;
arfidata = double(data.arfidata);
clear data;  
[x1 x2 x3] = ndgrid(xt,xt,xp);
dx = (x2-x1).*sign(0.5*(x1+x2)-x3);
DX = single(repmat(permute(dx,[4 1 2 3]),[size(arfidata,1),1,1,1]));
tic
if usejava('jvm')
    H = waitbar(0,sprintf('Cross-Correlating depth %0.0f/%0.0f',0,size(arfidata,1)));
    set(H,'Name',sprintf('%0.0f%%',0))
end
for zidx = 1:size(arfidata,1);
    if usejava('jvm')
    waitbar(zidx/size(arfidata,1),H,sprintf('Cross-Correlating depth %0.0f/%0.0f',zidx,size(arfidata,1)));
    set(H,'Name',sprintf('%0.0f%%',zidx/size(arfidata,1)*100))
    end
parfor pushidx = 1:numpushloc
    dxt = xt - xp(pushidx);
    disp = squeeze(arfidata(zidx,:,:,pushidx))';
    maxDelay = [max(abs(dxt))/dxdtRange(1)];
    tMax = maxDelay + kerLen;
    if tMax>max(t);
    tsrch = t(1):(t(2)-t(1)):tMax;
    srch = padarray(disp,[length(tsrch)-length(t) 0],'post');
    else
    tsrch = t;
    srch = disp;
    end
    [dxx dtt] = meshgrid(dxt,tsrch);
    msk = (sign(dxx).*dxx./dtt)<dxdtRange(2) & ((sign(dxx).*dxx./(dtt-kerLen))>dxdtRange(1) | dtt<kerLen);
    srch = srch.*msk;
    cc = nan(size(srch,1)*2-1,length(xt),length(xt));
    for refidx = 1:length(xt);
        ker = srch(:,refidx);
        tgtidx = find(sign(dxt(refidx))~=(-1*sign(dxt)) & abs(dxt(refidx))<abs(dxt));
        if any(tgtidx)
        cc(:,refidx,tgtidx) = xcorr2(srch(:,tgtidx),ker);
        end
    end
    lags = [-(size(srch,1)-1):(size(srch,1)-1)]*(t(2)-t(1));
    [CC(zidx,:,:,pushidx) DT(zidx,:,:,pushidx) pkMtx] = subsamplepeak(lags,cc,1);
    if saveQ
        Q(zidx,:,:,pushidx,:) = pkMtx;
    end
end
end

if usejava('jvm')
    close(H)
end


    