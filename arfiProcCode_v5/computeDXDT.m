function [DXDT CC] = computeDXDT(DX,t,sweidata,sws);
if ~exist('sws','var')
    sws = struct(...
        'kernel',[2 4 8] ...
        ,'dxdt_search',[0.5:0.1:8] ...
        ,'method','poly'...
        ,'smurf',1 ...
        ,'xrange',[0.25 4] ...
        ,'trange',[0 2.5] ...
        );
end

if sws.smurf;
    sweidata = permute(sweidata,[1 4 3 2]);
end

if length(unique(diff(t)))>1;
    t0 = t(:);
    sweidata0 = reshape(permute(sweidata,[3 1 2 4]),length(t),[]);
    t = (t0(1):min(unique(diff(t))):t0(end))';
    sweidataup = interp1(t0,sweidata0,t);
    sweidata = permute(reshape(sweidataup,length(t),size(sweidata,1),size(sweidata,2),size(sweidata,4)),[2 3 1 4]);
end

[vel, t1] = differentiateDisplacements(sweidata,t,2000);
vel = removeLateralMedian(vel);
[vL, vR] = LeftRightFilter(vel);
DX3 = repmat(permute(DX,[1 2 4 3]),[1 1 size(vel,3) 1]);
vLR = vel;
vLR(DX3<0) = 2*vL(DX3<0);
vLR(DX3>0) = 2*vR(DX3>0);
T3 = repmat(reshape(t1(:),1,1,[]),[size(DX3,1),size(DX3,2),1,size(DX3,4)]);
vFilt = vLR;
vFilt(abs(DX3)<sws.xrange(1)) = nan;
vFilt(abs(DX3)>sws.xrange(2)) = nan;
vFilt(T3<sws.trange(1)) = nan;
vFilt(T3>sws.trange(2)) = nan;

numloc = size(vFilt,4);
DXDT = nan(size(vFilt,1),numloc-min(sws.kernel)+1,numloc,length(sws.kernel));
T0 = DXDT;
CC = DXDT;
for kidx = 1:length(sws.kernel)
    kernel = sws.kernel(kidx);
    switch sws.method
        case 'poly'
            [pk TTP0] = subsamplepeak(t1,vFilt,3);
            for sourceidx = 1:numloc
                TTP  = TTP0(:,:,1,sourceidx);
                clear TTPtmp DXtmp
                for i = 1:kernel
                    TTPtmp(:,:,i) = TTP(:,(0:size(TTP,2)-kernel)+i);
                    DXtmp(:,:,i) = abs(DX(:,(0:size(TTP,2)-kernel)+i,sourceidx));
                end
                DTDX = nansum((DXtmp-repmat(nanmean(DXtmp,3),[1 1 kernel])).*(TTPtmp-repmat(nanmean(TTPtmp,3),[1 1 kernel])),3)./nansum(~isnan(TTPtmp).*(DXtmp-repmat(nanmean(DXtmp,3),[1 1 kernel])).^2,3);
                T0(:,1:size(DTDX,2),sourceidx,kidx) = nanmean(TTPtmp,3)-DTDX.*nanmean(DXtmp,3);
                DXDT(:,1:size(DTDX,2),sourceidx,kidx) = 1./DTDX;
            end
        case 'cc0'
            fprintf('Computing shifts and cross correlation coefficients:\n')
            fprintf('0%%%s100%%\n  ',95*ones(1,numloc))
            fprintf('\n')
            lags = [0:30];
            dxdt = nan(size(DXDT,1),size(DXDT,2),size(DXDT,3));
            cc = nan(size(DXDT,1),size(DXDT,2),size(DXDT,3));
            parfor idx = 1:numloc;
                fprintf('\b>\n');
                DXi = DX(:,:,idx);
                disp = single(vFilt(:,:,:,idx));
                disp = (disp-repmat(nanmean(disp,1),[size(disp,1) 1 1]))./repmat(std(disp,1,1),[size(disp,1) 1 1]);
                dispnormref = disp(:,(1:numloc-kernel+1),:);
                dispnormtgt = disp(:,(kernel:numloc),[1 1:end end]);
                dispnormtgt(:,:,[1 end]) = nan;
                N = length(t1);
                tau = t1(2)-t1(1);
                lags = [-30:30];
                CCcoef = zeros(size(dispnormref,1),size(dispnormref,2),length(lags));
                for lagidx = 1:length(lags)
                    tidx = (1:N)+lags(lagidx)+1;
                    CCcoef(:,:,lagidx) = nanmean(dispnormref.*dispnormtgt(:,:,max(1,min(N+2,tidx))),3);
                end
                [ccpk, tpk] = subsamplepeak(lags,CCcoef,3);
                dt = nan(size(DXDT,1),size(DXDT,2));
                dt(:,1:numloc-kernel+1) = tau*squeeze(tpk);
                dx = nan(size(DXDT,1),size(DXDT,2));
                dx(:,1:numloc-kernel+1) = DXi(:,(1:numloc-kernel+1))-DXi(:,(kernel:numloc));
                sgn = nan(size(DXDT,1),size(DXDT,2));
                sgn(:,1:numloc-kernel+1) = -1*sign((DXi(:,(1:numloc-kernel+1))+DXi(:,(kernel:numloc)))/2);
                cctmp = nan(size(DXDT,1),size(DXDT,2));
                cctmp(:,1:numloc-kernel+1) = ccpk;
                cc(:,:,idx) = cctmp;
                dxdt(:,:,idx) = sgn.*dx./dt;
            end
            DXDT(:,:,:,kidx) = dxdt;
            CC(:,:,:,kidx) = cc;
            fprintf('\n')
        case 'cc'
            fprintf('Computing shifts and cross correlation coefficients:\n')
            fprintf('0%%%s100%%\n  ',95*ones(1,numloc))
            fprintf('\n')
            lags = [0:30];
            dxdt = nan(size(DXDT,1),numloc-kernel+1,size(DXDT,3));
            cc = nan(size(DXDT,1),numloc-kernel+1,size(DXDT,3));
            parfor idx = 1:numloc;
                fprintf('\b>\n');
                DXi = DX(:,:,idx);
                disp = single(vFilt(:,:,:,idx));
                disp = (disp-repmat(nanmean(disp,1),[size(disp,1) 1 1]))./repmat(std(disp,1,1),[size(disp,1) 1 1]);
                dispnormref = disp(:,(1:numloc-kernel+1),:);
                dispnormtgt = disp(:,(kernel:numloc),:);
                %dispnormtgt(:,:,[1 end]) = nan;
                N = length(t1);
                tau = t1(2)-t1(1);
                lags = [-30:30];
                %dx = nan(size(DXDT,1),size(DXDT,2));
                dx = DXi(:,(1:numloc-kernel+1))-DXi(:,(kernel:numloc));
                sgn = -1*sign((DXi(:,(1:numloc-kernel+1))+DXi(:,(kernel:numloc)))/2);
                CCcoef = zeros(size(dispnormref,1),size(dispnormref,2),length(lags));
                ct_test =  [0.5:0.1:8];
                for testidx = 1:length(ct_test);
                dti = dx./ct_test(testidx)/tau;
                dispnormtgtshift = shift3(dispnormtgt,sgn.*dti);
%                 subplot(131)
%                 imagesc(squeeze(dispnormtgtshift(40,:,:)),[0 1])
%                 subplot(132)
%                 imagesc(squeeze(dispnormref(40,:,:)),[0 1]);
%                 subplot(133)
%                 imagesc(squeeze(dispnormref(40,:,:).*dispnormtgtshift(40,:,:)),[0 1])
%                 pause(0.1);
                CCcoef(:,:,testidx) = nanmean(dispnormref.*dispnormtgtshift,3);
                end
                [ccpk, ctpk] = subsamplepeak(ct_test,CCcoef,3);
                cc(:,:,idx) = ccpk;
                dxdt(:,:,idx) = ctpk;
            end
            DXDT(:,1:size(dxdt,2),:,kidx) = dxdt;
            CC(:,1:size(cc,2),:,kidx) = cc;
            fprintf('\n')            
    end
    
    
    
end
