ecgdatadir = '/getlab/pjh7/SC2000/ECGData/';
filepaths = {'/getlab/pjh7/SC2000/data/20130530_RV_Mmode/','/getlab/pjh7/SC2000/data/20130530_RA_Mmode/','/getlab/pjh7/SC2000/data/20130530_LV_Mmode/'};

for filepathidx = 1:3
    filepath = filepaths{filepathidx};
%fileindices = 1:2:50;
SWIF_files = dir(fullfile(filepath,'SWIF_AData*.bin'));
fileindices = 1:length(SWIF_files);
method = 'PesaventoFlux';
tic
for fileidx = fileindices
clear M im sp
timestamp = SWIF_files(fileidx).name(regexp(SWIF_files(fileidx).name,repmat('\d',1,14))+[0:13]);
resfile = fullfile(filepath,'res',method,sprintf('res_all_%s.mat',timestamp));
if ~exist(resfile,'file')
    continue
end
par = load(fullfile(filepath,sprintf('par_%s.mat',timestamp)));
toc1 = toc;
fprintf('Loading %s...',resfile)
res = load(resfile);
fprintf('done (%0.2fs)\n',toc-toc1);
ecgfile = fullfile(filepath,'ecg',sprintf('ecg_%s.mat',timestamp));
if exist(ecgfile,'file')
    ecg = load(ecgfile);
end

toc1 = toc;
fprintf('Filtering...')
v0 = diff(double(res.arfidata0_all)*res.arfi_scale,1,3);
v0(:,:,56,1:end-1) = 0.5*(v0(:,:,55,1:end-1)+v0(:,:,1,2:end));
v0 = v0(:,:,:);
v0 = v0(:,:,[end 1:end-1]);
[v2 t1] = differentiateDisplacements(cumsum(v0,3),res.t,[100 1000]);
v2 = reshape(v2(:,:,[1:end end]),size(res.arfidata0_all));
v2 = v2(:,:,1:end-1,:);
clear v0
[vL vR] = LeftRightFilter(v2-repmat(median(v2,2),[1 32 1 1]));
clear v2
vLR = [vL(:,1:16,:,:) vR(:,17:32,:,:)];
%fprintf('%3.0f/%3.0f...',0,192)
% [v t1] = differentiateDisplacements(double(res.arfidata_all).*res.arfi_scale,res.t,[100 1000]);
% for i = 1:192
%     fprintf([11*ones(1,7) sprintf('%3.0f/%3.0f...',i,192)]);
% vel = v(:,:,:,i);
% vel = removeLateralMedian(vel);
% [vL vR] = LeftRightFilter(vel);
% vLR(:,:,:,i) = [vL(:,1:16,:) vR(:,17:32,:)];
% end
clear vL vR vel v
%vLR = reshape(vLR,size(vLR,1),32,[]);
fprintf('done (%0.2fs)\n',toc-toc1);
toc1 = toc;
[X Z DX] = sectorCoordinates(res.theta,res.theta,res.axial,res.apex);
kernel_z = 1;
xrange = [0 3];
kernel_length = 1/180; %s
search_region = [0.5 8]; %m/s
search_step = 0.05;
ct_search = search_region(1):search_step:search_region(2);
z_out = linspace(5,20);
t = (0:size(vLR,3)-1)*mean(diff(t1))*1e-3;
sws = zeros(length(z_out),length(t));
%fprintf('Depth %3.0f/%3.0f',0,length(z_out));
% fprintf('Normalizing %3.0f/%3.0f',0,192)
% for i = 1:192
%        fprintf([8*ones(1,7) sprintf('%3.0f/%3.0f',i,192)]);
% vLR(:,:,:,i) = vLR(:,:,:,i)-repmat(mean(vLR(:,:,:,i),3),[1 1 size(vLR,3) 1]);
% vLR(:,:,:,i) = vLR(:,:,:,i)-repmat(std(vLR(:,:,:,i),0,3),[1 1 size(vLR,3) 1]);
% end


k = repmat(permute(gausswin(7,1),[2 3 4 1]),[5 1 1 1]);
bmsk = repmat(permute((double(res.bmodedata_all)/255).^2,[1 2 4 3]),[1 1 size(vLR,3) 1]);
%vLR = convn(bmsk.*vLR,k,'same')./convn(bmsk,k,'same');
DX3 = repmat(DX,[1 1 length(t1)]);
T3 = repmat(permute(t1(:),[2 3 1]),[size(DX,1) size(DX,2) 1]);
msk = ((abs(DX3)./(T3-1.5))>1)|T3<1.5;
vLR = vLR.*repmat(msk,[1 1 1 size(vLR,4)]);
%%
% dxdt = [];
% fprintf('z: %3.0f/%3.0f t: %3.0f/%3.0f\n',0,size(vLR,1),0,size(vLR,4));
% ctest = 0.5:0.1:15;
% ncc = zeros(size(vLR,1),size(vLR,4),length(ctest));
% parfor zidx = 1:size(vLR,1) 
%     fprintf([8*ones(1,19) '%3.0f/%3.0f t: %3.0f/%3.0f\n'],zidx,size(vLR,1),0,size(vLR,4));
%     x = nanmean(DX(zidx,:),1);
%     ncc0 = zeros(size(vLR,4),length(ctest));
%     %for tidx = 1:size(vLR,4)
%         fprintf([8*ones(1,8) '%3.0f/%3.0f\n'],tidx,size(vLR,4));
%         slice = squeeze(vLR(zidx,:,:,:));
%         v = squeeze(mean(slice(floor((size(slice,1)+1)/2):ceil((size(slice,1)+1)/2),:,:),1));
%         tst = zeros(length(x),length(t1),size(vLR,4),length(ctest));
%         for i = 1:length(ctest)
%             dt = (abs(x)/ctest(i));
%             for j = 1:length(x);
%                 tst(j,:,:,i) = interp1(t1,v,t1-dt(j));
%             end
%         end
%         ncc00 = zeros(1,size(vLR,4),length(ctest));
%         for i = 1:length(ctest)
%             ncc0(1,:,i) = nanmean(nanmean(tst(:,:,:,i).*slice,1),2);
%         end
%     ncc(zidx,:,:) = ncc0;
% end
% [pk sws] = subsamplepeak(ctest,ncc,3);
% fprintf('\n');

%%

refidx = [1:15 18:32];
ccmethod = 'progressive';
switch ccmethod
    case 'progressive'
        tgtidx = [2:31];
        shifts = [0:10];
    case 'anchored'        

        tgtidx = [16*ones(1,15) 17*ones(1,15)];
        shifts = [0:20];

end

vTgt = vLR(:,tgtidx,:,:);
vRef = vLR(:,refidx,:,:);
fprintf('Computing cross-correlation coefficients...')
fprintf('%2.0f/%2.0f',0,length(shifts))
CC = nan(size(vLR,1),size(vTgt,2),length(shifts),size(vLR,4));
for shiftidx = 1:length(shifts)
    fprintf([8*ones(1,5) sprintf('%2.0f/%2.0f',shiftidx,length(shifts))]);
    shift = shifts(shiftidx);
    %vTgtShift = reshape(vTgt(:,:,[ones(1,shift) 1:size(vTgt,3)-shift]),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vRefShift = reshape(vRef(:,:,[shift+1:end size(vRef,3)*ones(1,shift)]),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vTgtShift = reshape(vTgt(:,:,max(1,min(size(vTgt,3),(1:size(vTgt,3))-shift))),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vRefShift = reshape(vRef(:,:,max(1,min(size(vTgt,3),(1:size(vRef,3))+shift))),size(vTgt,1),size(vTgt,2),length(t1),[]);
    vTgtShift = vTgt(:,:,max(1,min(size(vTgt,3),(1:size(vTgt,3))-shift)),:);
    vRefShift = vRef(:,:,max(1,min(size(vTgt,3),(1:size(vRef,3))+shift)),:);
    
    %vTgtShift = vTgtShift-repmat(mean(vTgtShift,3),[1 1 size(vTgtShift,3) 1]);
    %vRefShift = vRefShift-repmat(mean(vRefShift,3),[1 1 size(vRefShift,3) 1]);
    %vTgtShiftSig = std(vTgtShift,1,3);
    %vRefShiftSig = std(vRefShift,1,3);
    %imagesc([squeeze(vTgtShift(70,:,:,1));squeeze(vRefShift(70,:,:,1))],[-5 5])
    %pause(0.1)
    CC(:,:,shiftidx,:) = mean(vTgtShift.*vRefShift,3);%./(size(vTgtShift,3).*vTgtShiftSig.*vRefShiftSig);
    %plot(shifts,squeeze(CC(70,:,:,2))','.-')
end
fprintf('...done\n')
[cc_coef dt] = subsamplepeak(mean(diff(t1))*shifts,CC,3);
cc_coef = squeeze(cc_coef);
dx = (DX(:,refidx)-DX(:,tgtidx)).*sign(DX(:,refidx));
x = (DX(:,refidx)+DX(:,tgtidx))/2;
dxdt = squeeze(repmat(dx,[1 1 1 size(dt,4)])./dt);
%%
% clf
% subplot(121)
% im = imsurf(res.theta,res.axial,res.apex,vLR(:,:,1,1),[-5 5]);
% hold on
% p(1) = plot(t1(1)*squeeze(nanmedian(dxdt(:,:,1),2)),res.axial,'g-');
% p(2) = plot(t1(1)*squeeze(nanmedian(-1*dxdt(:,:,1),2)),res.axial,'g-');
% subplot(122)
% im1 = imagesc([0:191]*(1/180),res.axial,squeeze(nanmedian(dxdt,2)),[-8 8]);
% set(im1,'alphaData',max(0,min(1,3*squeeze(nanmedian(bimg(:,14:19,:).^2,2)))))
% hold on
% p(3) = plot([0 0],res.axial([1 end]),'g-');
% for j = 1:100;
%     set(p(3),'xdata',(j-1)*(1/180)*[1 1]);
%     for i = 5:30;
%         tic;
%         set(im,'cdata',padarray(vLR(:,:,i,j),[1 1],'post'));
%         set(p(2),'xdata',t1(i)*squeeze(nanmedian(dxdt(:,:,j),2)));
%         set(p(1),'xdata',t1(i)*squeeze(nanmedian(-1*dxdt(:,:,j),2)));
%         drawnow;
%         while(toc<0.02);
%         end;
%     end;
% end


%%
dxdtpath = fullfile(filepath,'dxdt',ccmethod);
if ~exist(dxdtpath,'dir')
    mkdir(dxdtpath);
end
dxdtfile = fullfile(dxdtpath,sprintf('dxdt_%s.mat',timestamp));
save(dxdtfile,'dxdt','cc_coef');


% 
% 
% for frameidx = 1:192
% for depthidx = 1:length(z_out)
% fprintf([8*ones(1,7) sprintf('%3.0f/%3.0f',depthidx,length(z_out))]);
% msk = repmat((Z>(z_out(depthidx)-kernel_z/2))&(Z<=(z_out(depthidx)+kernel_z/2)),[1 1 size(vLR,3)]);
% xx = squeeze(sum(msk(:,:,1).*X,1)./sum(msk(:,:,1),1))*1e-3;
% slice = squeeze(sum(msk.*vLR(:,:,:,frameidx),1)./sum(msk,1));
% 
% end
% end
% dslice = repmat(,[1 1 length(ct_search)]);
%     for j = 1:length(ct_search);
%         shift = abs(xx)/ct_search(j);
%         for k = 1:length(xx);
%             dslice(k,:,j) = interp1(t,dslice(k,:,j),t+shift(k),'linear',0);
%         end
%     end
%     dslice = dslice(:,:,:).*dslice([16*ones(1,16) 17*ones(1,16)],:,:);
%     xidx = abs(xx*1e3)>=xrange(1) & abs(xx*1e3)<=xrange(2);
%     [pk idx] = max(filter2(ones(56,1),squeeze(sum(dslice(xidx,:,:),1)),'same'),[],2);   
%     sws(depthidx,:) = ct_search(idx);
% end
end
end
