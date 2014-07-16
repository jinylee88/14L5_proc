function [dttps_file] = process_DTTPS(datapath,useparfor)
%process_DTTPS processes "complete" dataset

%% Inital Paths and parameters
numSWEI = 1;
numARFI = 0;
%%
tic            
addpath /luscinia/sjr6/SC2000/arfiProcCode/
addpath /raidScratch/abp19/14L5HarmonicFocusedUnfocusedTrack/arfiProcCode_Filters/
swifFiles = dir(fullfile(datapath,'SWIF_AData*.bin'));
RESFILE= dir(fullfile(datapath,'resS_*.mat'));
dimFiles = dir(fullfile(datapath,'SWIF_ADataDims*.txt'));

for ind = 1:length(RESFILE)
        timeStamp = RESFILE(ind).name(6:end-4);
        resfile = ['resS_' timeStamp '.mat'];
        parfile = ['parS_' timeStamp '.mat'];
        if ~exist(parfile,'file')
        fprintf('par file %s not found!\n',fullfile(pwd,parfile));
%        keyboard
        movefile(files(ind).name,['f' swifFiles(ind).name]);
        continue
        end
        if ~exist(resfile, 'file')
            par = load(parfile);
            par.nref = 1;
            par.ntrack = [39 1];
            par.ensemble = par.nref+par.npush*length(par.pushFocalDepth)+sum(par.ntrack);
            save(parfile, '-struct', 'par')
            procArfi(swifFiles(ind).name);
        else
            fprintf(1, '%s\n', resfile);
        end
end           
resFiles = dir(fullfile(datapath,'resS_*.mat'));
resFiltFiles = dir(fullfile(datapath,'resFilt_*.mat'));
parFiles = dir(fullfile(datapath,'parS_*.mat'));
dttps_file = fullfile(datapath,'dttps.mat');
ttpData_file = fullfile(datapath,'ttpData.mat');
%% Get SWIF data names

bData = [];
par = load(fullfile(datapath,parFiles(1).name));
nBeams = par.nBeams;
numBeamGroups = par.numBeamGroups;
bidx = [1:par.nref par.nref+par.npush*length(par.pushFocalDepth)+(1:par.ntrack(1))];
%% Load SWIF data
clear bmPattern
for resIdx = 1:numSWEI
    %[data swifParams]= readSwif(fullfile(datapath,swifFiles(resIdx).name),fullfile(datapath,dimFiles(resIdx).name));
    load(fullfile(datapath,parFiles(resIdx).name),'trackParams');
    bmPattern(resIdx,:) = trackParams.rxMultibeamParams.beamPatternP;
  %  bDatai = reshape(mean(db(abs(complex(data.I(:,:,bidx),data.Q(:,:,bidx)))),3),size(data.I,1),par.nBeams,[]);
  %  bData = cat(3,bData,bDatai);
end
clear bDatai data
%% Determine receive beam positions
bmPatternFull = zeros(size(bmPattern,1)*numBeamGroups,nBeams);
for resIdx = 1:numSWEI
    for pushIdx = 1:numBeamGroups
        bmPatternFull((resIdx-1)*numBeamGroups+pushIdx,:) = bmPattern(resIdx,:)+(pushIdx-1)+mod(resIdx-1,numSWEI/2)*par.numBeamGroups;
    end
end
bmIdxFull = bmPatternFull-min(bmPatternFull(:))+1;
%% Retrieve B-mode data
nPush = numBeamGroups*numSWEI;
szX = max(bmPatternFull(:))-min(bmPatternFull(:))+1;
% bData1 = nan(size(bData,1),szX,nPush/2);
% for i = 1:(nPush/2)
%     bData1(:,bmIdxFull(i+(nPush/2),:),i) = bData(:,:,i+(nPush/2));
%     bData1(:,bmIdxFull(i,:),i) = bData(:,:,i);
% end

%% Load SWI displacement data
clear res
for resIdx = 1:numSWEI;
    resFile = fullfile(datapath,resFiles(resIdx).name);
    fprintf('loading %s...',resFile)
    res(resIdx) = load(resFile);
    fprintf('\n');
end
%% Load ARFI displacement data
clear res_arfi
%for resIdx = (1:numARFI)
%    resFile = fullfile(datapath,resFiles(resIdx+numSWEI).name);
 %   fprintf('loading %s...',resFile)
 %   res_arfi(resIdx) = load(resFile);
 %   fprintf('\n');
%end
%% Combine data and get spatial dimensions
sweidata = [res.arfidata];
latMat = [res.lat];
axial = res(1).axial;
t = res(1).t;
dLatR = res(1).lat-res(1).lat(ones(1,nBeams),:);
%dLatL = res(end).lat-res(end).lat(nBeams*ones(1,nBeams),:);
clear res
%%
%normData = load(fullfile(datapath,'..','normDataDisp.mat'));
normData = load('/getlab/pjh7/SC2000/FullySampled12L4/data/20140108/3kPa/1p5mm/normDataDisp.mat');
%keyboard

%downsampleFactor = mean(diff(axial))/mean(diff(res_arfi.axial));
%res_arfi.normData = single(resample(double(normData.normData),1,round(downsampleFactor)));
%res_arfi.arfidata = res_arfi.arfidata(:,2:3:end,:);
%for i = 1:size(res_arfi.arfidata,3);
%    arfidata(:,:,i) = single(resample(double(res_arfi.arfidata(:,:,i)),1,round(downsampleFactor)));
%end
%res_arfi.arfidata = arfidata;
%res_arfi.axial = resample(double(res_arfi.axial),1,round(downsampleFactor));
%zaxidx = find((res_arfi.axial>3.75 & res_arfi.axial<9.25) | (res_arfi.axial>32));
%tidx = find(res_arfi.t<-0.25 | res_arfi.t>2.6);
%[resid axmotion] = linearmotionfilter(permute(res_arfi.arfidata,[3 2 1]),res_arfi.axial,zaxidx,1);
%axmotion = permute(axmotion,[3 2 1]);
%[resid motion] = linearmotionfilter(axmotion,res_arfi.t,tidx,1);
%res_arfi.arfidata_filt = res_arfi.arfidata-motion;
%normdata = repmat(permute(res_arfi.normData,[1 3 2]),[1 size(res_arfi.arfidata,2) 1]);
%[resid axmotion] = linearmotionfilter(permute(normdata,[3 2 1]),res_arfi.axial,zaxidx,1);
%axmotion = permute(axmotion,[3 2 1]);
%[resid motion] = linearmotionfilter(axmotion,res_arfi.t,tidx,1);
%normdata_filt = normdata-motion;
%res_arfi.arfidata_norm = res_arfi.arfidata_filt./normdata_filt;
%res_arfi.arfidata_norm(res_arfi.arfidata_norm<0) = 0;
%res_arfi.arfidata_norm_inv = 1./res_arfi.arfidata_norm;

%arfidata = res_arfi.arfidata_norm;

%arfidata = res_arfi.arfidata./repmat(permute(normData.normData,[1 3 2]),[1 size(res_arfi.arfidata,2), 1]);
%downsampleFactor = mean(diff(axial))/mean(diff(res_arfi.axial));
%for i = 1:size(arfidata,3);
%    arfidata_ds(:,:,i) = single(resample(double(arfidata(:,:,i)),1,round(downsampleFactor)));
%end
%axial_arfi = res_arfi.axial;%resample(double(res_arfi.axial),1,round(downsampleFactor));
%lat_arfi = res_arfi.lat(2:3:end);
%arfidata = arfidata_ds(:,2:3:end,:);
%clear arfidata_ds;
%% Differentiate displacements
toc1 = toc;
fprintf('Differentiating...');
%v = zeros(size(sweidata,1),size(sweidata,2),size(sweidata,3)-20,'single');
arfidata_tmp = sweidata(:,:,20:end);
t_tmp = t(20:end);
%[arfidata_tmp,t_tmp]=UpsampleTimeAndDispPlanes(arfidata_tmp,t_tmp,48);
v = zeros(size(sweidata,1),size(sweidata,2),length(t_tmp)-1,'single');
t1 = 0.5*t_tmp(1:end-1) + 0.5*t_tmp(2:end);
%%UPSAMPLE HERE TO 50KHz

if useparfor
parfor i = 1:size(sweidata,2)
    v(:,i,:) = differentiateDisplacements(arfidata_tmp(:,i,:),t_tmp,[50 1000]);
end    
else 
for i = 1:size(sweidata,2)
    v(:,i,:) = differentiateDisplacements(arfidata_tmp(:,i,:),t_tmp,[50 1000]);
end
end
clear arfidata_tmp t_tmp sweidata
fprintf('done (%0.2fs)\n',toc-toc1);
%% Create left-right masks
[Z DXR DT1] = ndgrid(axial,dLatR(:,1),t1);
%[Z DXL DT1] = ndgrid(axial,dLatL(:,1),t1);
DZ = Z-par.pushFocalDepth(1);
MSKR = (abs(DXR+1)./(DT1))>0.5 & (abs(DXR)./(DT1+.05))<6 & (0.8*abs(DZ)-4*DXR)<1;
%MSKL = (abs(DXL-1)./(DT1))>0.5 & (abs(DXL)./(DT1+.05))<6 & (0.8*abs(DZ)+4*DXL)<1;
MSKR = (abs(DXR+1)./(DT1))>0.5 & (abs(DXR)./(DT1+.05))<6 ;
%MSKL = (abs(DXL-1)./(DT1))>0.5 & (abs(DXL)./(DT1+.05))<6 ;

%% Left-Right filter
vLR = v;
toc1 = toc;
fprintf('Left-Right Filtering...');
%H = waitbar(0,'Left-Right Filtering...');
for i = 1:nPush;
%    waitbar(i/nPush,H);
    [vL vR] = LeftRightFilter(v(:,(i-1)*nBeams+(1:nBeams),:));
    vR(~MSKR) = nan;
    vLR(:,(i-1)*nBeams+(1:nBeams),:) = vR;
end
% for i = nPush/2 + [1:nPush/2];
% %    waitbar(i/nPush,H);
%     [vL vR] = LeftRightFilter(v(:,(i-1)*nBeams+(1:nBeams),:));
%     vL(~MSKL) = nan;
%     vLR(:,(i-1)*nBeams+(1:nBeams),:) = vL;
% end
fprintf('done (%0.2fs)\n',toc-toc1);
%clear v vL vR
%% Left-Right filter (Push)
vLR1 = v*0;
%% 
Idx = zeros(1,size(v,2));
for TrackIdx = 1:(nPush + nBeams-1);
    idx = ((TrackIdx-1)*nBeams)+1:-1*(nBeams-1):TrackIdx;
    if length(idx)>nBeams
        idx = idx(1:nBeams);
    end
    %idx = idx(idx<=size(v,2)/2);
    idx = idx(idx<=size(v,2));
    Idx(idx) = Idx(idx)+1;
    [vL vR] = LeftRightFilter(v(:,idx,:));
    vLR1(:,idx,:) = vR;
end
for PushIdx = 1:nPush
    vtmp = vLR1(:,(PushIdx-1)*nBeams+[1:nBeams],:);
    vtmp(~MSKR) = nan;
    vLR1(:,(PushIdx-1)*nBeams+[1:nBeams],:) = vtmp;
end
% for TrackIdx = (1:nPush/2 + nBeams-1);
%     idx = size(v,2)/2+((TrackIdx-1)*nBeams)+1:-1*(nBeams-1):(size(v,2)/2);
%     if length(idx)>nBeams
%         idx = idx(1:nBeams);
%     end
%     idx = idx(idx>size(v,2)/2 & idx<=size(v,2));
%     Idx(idx) = Idx(idx)+1;
%     [vL vR] = LeftRightFilter(v(:,idx,:));
%     vLR1(:,idx,:) = vL;
% end
% for PushIdx = (nPush/2)+(1:nPush/2)
%     vtmp = vLR1(:,(PushIdx-1)*nBeams+[1:nBeams],:);
%     vtmp(~MSKL) = nan;
%     vLR1(:,(PushIdx-1)*nBeams+[1:nBeams],:) = vtmp;
% end
clear v vtmp vL vR
%close(H);
%% Axially median filter
toc1 = toc;
fprintf('Median Filtering...')
if useparfor
    parfor i = 1:size(vLR,2) %parfor
        vLR(:,i,:) = medfilt1(double(vLR(:,i,:)),9,[],1);
        vLR1(:,i,:) = medfilt1(double(vLR1(:,i,:)),9,[],1);
    end
    else
    for i = 1:size(vLR,2) %parfor
        vLR(:,i,:) = medfilt1(double(vLR(:,i,:)),9,[],1);
        vLR1(:,i,:) = medfilt1(double(vLR1(:,i,:)),9,[],1);
    end
end
% if useparfor
%     parfor i = 1:size(arfidata,2) %parfor
%         arfidata(:,i,:) = medfilt1(double(arfidata(:,i,:)),9,[],1);
%     end
% else
%     for i = 1:size(arfidata,2) %parfor
%         arfidata(:,i,:) = medfilt1(double(arfidata(:,i,:)),9,[],1);
%     end
% end
fprintf('done (%0.2fs)\n',toc-toc1);

%% Find time-to-peak slope
toc1 = toc;
fprintf('Finding Peaks...');
[pk pk1 tpk tpk1] = deal(zeros(size(vLR,1),size(vLR,2)));
if useparfor
    parfor i = 1:size(vLR,2) %parfor
        [pk(:,i) tpk(:,i)] = subsamplepeak(t1,vLR(:,i,:),3);
        [pk1(:,i) tpk1(:,i)] = subsamplepeak(t1,vLR1(:,i,:),3);
    end
else
    for i = 1:size(vLR,2) %parfor
        [pk(:,i) tpk(:,i)] = subsamplepeak(t1,vLR(:,i,:),3);
        [pk1(:,i) tpk1(:,i)] = subsamplepeak(t1,vLR1(:,i,:),3);
    end
end
fprintf('done (%0.2fs)\n',toc-toc1);
tpk = reshape(tpk,size(vLR,1),nBeams,[]);
tpk1 = reshape(tpk1,size(vLR,1),nBeams,[]);

clear vLR vLR1
%% Assign ttps data to single matrix
[tpkL tpkL1 tpkR tpkR1] = deal(nan(size(tpk,1),szX,nPush/2));
dx = median(diff(latMat(2:end,1)));
%dxL = nan(1,szX,nPush/2);
dxR = nan(1,szX,nPush);
for i = 1:nPush
%    tpkL(:,bmIdxFull(i+(nPush/2),:),i) = tpk(:,:,i+(nPush/2));
%    tpkL1(:,bmIdxFull(i+(nPush/2),:),i) = tpk1(:,:,i+(nPush/2));
%    dxL(1,bmIdxFull(i+(nPush/2),:),i) = dx*bmPattern(end,:);
    tpkR(:,bmIdxFull(i,:),i) = -1*tpk(:,:,i);
    tpkR1(:,bmIdxFull(i,:),i) = -1*tpk1(:,:,i);
    dxR(1,bmIdxFull(i,:),i) = dx*bmPattern(1,:);
end
tpkLR = tpkR;
tpkLR1 = tpkR1;

xr = [0:size(tpkLR,2)-1]*dx;xr = xr-mean(xr);
xs = [0:size(tpkLR,3)-1]*dx;xs = xs-mean(xs);
%[Z XR XS] = ndgrid(axial,xr,xs);
clear tpk tpk1
%% Remove offset and combine estimates from overlapping pushes/tracks
% zKrnl = 3;
% xKrnl = 9;
% dtpkLRstl = diff(tpkLR,1,3);
% dttp_stl = zeros(size(dtpkLRstl,1),size(dtpkLRstl,3));
% parfor i = 1:size(dtpkLRstl,1)
%     rowmed = nan(1,size(dtpkLRstl,3));
%     for j = 1:size(dtpkLRstl,3);
%         ttp_tmp = dtpkLRstl(max(1,min(size(dtpkLRstl,1),i+[-(zKrnl-1)/2:(zKrnl-1)/2])),:,max(1,min(size(dtpkLRstl,3),j+[-(xKrnl-1)/2:(xKrnl-1)/2])));
%         rowmed(j) = nanmedian(ttp_tmp(:));
%     end
%     dttp_stl(i,:) = rowmed;
% end
%
% offsets = -1*nanmedian(diff(tpkLR,1,3),2);
% tpkLR1 = tpkLR;
% for i = 2:size(tpkLR1,3);
%     tpkLR1(:,:,i) = tpkLR(:,:,i)+repmat(sum(offsets(:,:,1:i-1),3),[1 size(tpkLR,2)]);
% end
%
% dtpkLRmtl = -1*diff(tpkLR1,1,2);
% dttp_mtl = zeros(size(dtpkLRmtl,1),size(dtpkLRmtl,2));
% parfor i = 1:size(dtpkLRmtl,1)
%     rowmed = nan(1,size(dtpkLRmtl,2));
%     for j = 1:size(dtpkLRmtl,2);
%         ttp_tmp = dtpkLRmtl(max(1,min(size(dtpkLRmtl,1),i+[-(zKrnl-1)/2:(zKrnl-1)/2])),max(1,min(size(dtpkLRmtl,2),j+[-(xKrnl-1)/2:(xKrnl-1)/2])),:);
%         rowmed(j) = nanmedian(ttp_tmp(:));
%     end
%     dttp_mtl(i,:) = rowmed;
% end

%% Clip data so that MTL and STL have same spatial sampling
% tpkRclip = tpkR(:,(bmIdxFull(1,1):bmIdxFull(nPush,1))-bmPattern(1,1),:);
% tpkLclip = tpkL(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
% tpkR1clip = tpkR1(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
% tpkL1clip = tpkL1(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
% dxLclip =  dxL(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
%dxRclip =  dxR(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
dlat =  dxR;
%tpkLRclip = tpkLR(:,(bmIdxFull(1,1):bmIdxFull(nPush/2,1))-bmPattern(1,1),:);
tpkRclip=tpkR;
tpkR1clip=tpkR1;
%% "Align" ttps data
%offsetsL = -1*nanmean(diff(tpkLclip,1,3),2);
% tpkL1 = tpkLclip;
% for i = 2:size(tpkL1,3);
%     tpkL1(:,:,i) = tpkLclip(:,:,i)+repmat(nansum(offsetsL(:,:,1:i-1),3),[1 size(tpkLclip,2)]);
% end
%clear offsetsL
%tpkL1 = tpkL1-repmat(tpkL1(:,bmIdxFull(nPush/2),end),[1 size(tpkL1,2) size(tpkL1,3)]);
%ttpsL = -1*squeeze(nanmedian(tpkL1,3));
offsetsR = -1*nanmean(diff(tpkRclip,1,3),2);
tpkR1 = tpkRclip;
for i = 2:size(tpkR1,3);
    tpkR1(:,:,i) = tpkRclip(:,:,i)+repmat(nansum(offsetsR(:,:,1:i-1),3),[1 size(tpkRclip,2)]);
end
clear offsetsR
%tpkR1 = tpkR1-repmat(tpkR1(:,bmIdxFull(1,1),1),[1 size(tpkR1,2) size(tpkR1,3)]);
%ttpsL = -1*squeeze(nanmedian(tpkL1,3));
ttpsR = -1*squeeze(nanmedian(tpkR1,3));
%%
ttps =tpkRclip;
ttps1 = tpkR1clip;
%ttps1=ttps;

%% for now Adam Addded a temporary section while I figure out how to filter the data better. This will Nan any value in the ttps matrix that is euqal to t1(1) because this is caused by some filtering problem
ttps(ttps==t1(1) | ttps==-t1(2))=nan;
ttps1(ttps1==t1(1) | ttps==-t1(2))=nan;

ttps(ttps==-t1(1) | ttps==-t1(2))=nan;
ttps1(ttps1==-t1(1) | ttps==-t1(2))=nan;

clear tpkL1 tpkR1 tpkLR tpkLclip tpkRclip

%% Combine ttps data from overlapping pushes/tracks
dttpsRmtl = diff(ttpsR,1,2);
%dttpsLmtl = diff(ttpsL,1,2);
dttpsRstl = squeeze(nanmedian(diff(tpkR,1,3),2));
%dttpsLstl = squeeze(nanmedian(diff(tpkL,1,3),2));

%%
x1 = [0:nPush-1]*dx;
x1 = x1-mean(x1);
z1 = axial;
x=x1;
z=z1;
%x = lat_arfi;
%z = axial_arfi;
ttps(ttps==0)=nan;
ttps1(ttps1==0)=nan;
fprintf('saving %s...',dttps_file)
save(dttps_file,'ttps','ttps1','dlat','dttpsRmtl','dttpsRstl','x1','z1','x','z')
fprintf('done\n');

