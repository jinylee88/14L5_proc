function dxdtFile = progressiveNCC(resFileName,mode);
smurf_flag = ~isempty(strfind(mode,'smurf'));
if ~isempty(strfind(mode,'anchored'))
    ccmethod = 'anchored';
else
    ccmethod = 'progressive';
end
tic
filepath = fileparts(fileparts(fileparts(resFileName)));
if isempty(filepath)
    filepath = '../../';
end
timestamp = RetrieveTimeStamp(resFileName);
par = load(fullfile(filepath,sprintf('par_%s.mat',timestamp)));

toc1 = toc;
fprintf('Loading %s...',resFileName)
res = load(resFileName);
if ~isfield(res,'vTest')
    dxdtFile = '';
    fprintf('no vTest. Closing...\n');
    return
end
fprintf('done (%0.2fs)\n',toc-toc1);
toc1 = toc;
[pushbeamSort jj] = sort(par.pushbeamNum);
res.vTest = res.vTest(:,:,:,jj);
if smurf_flag
res.vTest = permute(res.vTest,[1 4 3 2]);
end
t1 = filter2([0.5 0.5],res.t,'valid');
for i = 1:length(par.pushbeamNum);
    res.dtheta(i,:) = res.theta-res.theta(i);
end
[X Z DX] = sectorCoordinates(res.theta,res.dtheta,res.axial,res.apex);
kernel_z = 1;
xrange = [0 3];
kernel_length = 1/180; %s
search_region = [0.5 8]; %m/s
search_step = 0.05;
ct_search = search_region(1):search_step:search_region(2);
k = repmat(permute(gausswin(7,1),[2 3 4 1]),[5 1 1 1]);
bmsk = repmat((double(res.bmodedata0)/255).^2,[1 1 size(res.vTest,3)]);
T3 = repmat(permute(t1(:),[2 3 1]),[size(DX,1) size(DX,2) 1]);
theta1 = filter2([0.5 0.5],res.theta,'valid');
dxdt = zeros(size(res.vTest,1),length(theta1),size(res.vTest,4));
cc_coef = dxdt;

fprintf('Computing cross-correlation coefficients...')
fprintf('%2.0f/%2.0f',0,size(res.vTest,4))
for pushIdx = 1:size(res.vTest,4);
fprintf([8*ones(1,5) sprintf('%2.0f/%2.0f',pushIdx,size(res.vTest,4))]);

    %fprintf('Push %0.0f/%0.0f...\n',pushIdx,size(res.vTest,4));
DX3 = repmat(DX(:,:,pushIdx),[1 1 length(t1)]);
msk = (((abs(DX3)./(T3-1.5))>1)|T3<1.5) & abs(DX3)<3;
vLR = double(res.vTest(:,:,:,pushIdx))*res.arfi_scale.*msk;
vLR(abs(DX3)>3) = nan;
refidx = [1:pushIdx-1 pushIdx+1:size(vLR,2)];
switch ccmethod
    case 'progressive'
        tgtidx = [2:pushIdx pushIdx:size(vLR,2)-1];
        shifts = [-5:10];
        vTgt = vLR(:,tgtidx,:);
        vRef = vLR(:,refidx,:);
    case 'anchored'        
        tgtidx = pushIdx*ones(1,size(vLR,2)-1);
        shifts = [-2:20];
        vTgt = repmat(nanmedian(vLR(:,max(1,min(size(vLR,2),pushIdx+[-1:1])),:),2),[1 size(vLR,2)-1 1]);
        vRef = vLR(:,refidx,:);

end
CC = nan(size(vLR,1),size(vTgt,2),length(shifts),size(vLR,4));
for shiftidx = 1:length(shifts)
    shift = shifts(shiftidx);
    %vTgtShift = reshape(vTgt(:,:,[ones(1,shift) 1:size(vTgt,3)-shift]),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vRefShift = reshape(vRef(:,:,[shift+1:end size(vRef,3)*ones(1,shift)]),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vTgtShift = reshape(vTgt(:,:,max(1,min(size(vTgt,3),(1:size(vTgt,3))-shift))),size(vTgt,1),size(vTgt,2),length(t1),[]);
    %vRefShift = reshape(vRef(:,:,max(1,min(size(vTgt,3),(1:size(vRef,3))+shift))),size(vTgt,1),size(vTgt,2),length(t1),[]);
    vTgtShift = vTgt(:,:,max(1,min(size(vTgt,3),(1:size(vTgt,3)))),:);
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
%fprintf('...done\n')
[cc_coef0 dt] = subsamplepeak(mean(diff(t1))*shifts,CC,3);
cc_coef(:,:,pushIdx) = squeeze(cc_coef0);
dx = (DX(:,refidx,pushIdx)-DX(:,tgtidx,pushIdx)).*sign(DX(:,refidx,pushIdx));
x = (DX(:,refidx,pushIdx)+DX(:,tgtidx,pushIdx))/2;
dt(cc_coef0<=0) = nan;
dt(dt<=0) = nan;
dt = stdFilt(dt,[11 3],1);
dxdt(:,:,pushIdx) = dx./dt;
axial = res.axial;
theta = res.theta;
apex = res.apex;
%%
% vTest(:,:,:,pushIdx) = vLR;
% figure(2);clf
% imagesc(t1,DX(66,:,pushIdx),reshape(median(double(vTest(66,:,:,pushIdx)),1),32,[]),[-20 20]);
% hold on
% p = plot(dt(66,:),dx(66,:).*sign(DX(66,refidx,pushIdx)),'.-');

end
fprintf('...done\n')
dxdtpath = fullfile(filepath,'dxdt',mode);
if ~exist(dxdtpath,'dir')
    mkdir(dxdtpath);
end
dxdtFile = fullfile(dxdtpath,sprintf('dxdt_%s.mat',timestamp));
save(dxdtFile,'axial','theta','theta1','DX','dxdt','cc_coef','apex');
end

