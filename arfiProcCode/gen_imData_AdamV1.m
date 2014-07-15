function imData = gen_imData_AdamV1(acqPath,regressionFilters,medianFilters,DSfactor, TrackLines, D2filt)
TrackBeams=TrackLines;
if ~exist('maxDX','var')
    maxDX = inf;
end
modalities = {'SWEI_MTL','SWEI_STL'}
%dttps_file = fullfile(acqPath,['dttpsTrack' num2str(TrackBeams) 'DS' num2str(DSfactor) '.mat']);
dttps_file = fullfile(acqPath,['dttps.mat']);
if ~exist(dttps_file,'file');
    [success output] = system('hostname');
    process_DTTPS(acqPath,strcmpi(deblank(output),'gudenaa.egr.duke.edu'));
end


clear bmodedata
fprintf('loading %s...',dttps_file);
load(dttps_file);
x = x-(mean(x));
fprintf('done\n');
dx = mean(diff(x));

%Axially average TTPS data

if ~exist('bmodedata','var')
    %% Get SWIF data names
    swifFiles = dir(fullfile(acqPath,'SWIF_AData*.bin'));
    dimFiles = dir(fullfile(acqPath,'SWIF_ADataDims*.txt'));
    % parFiles = dir(fullfile(acqPath,['parTrack' num2str(TrackBeams) 'DS' num2str(DSfactor) '_*.mat']));
    parFiles = dir(fullfile(acqPath,['parS_*.mat']));
    bData = [];
    par = load(fullfile(acqPath,parFiles(end).name),'c','fs');
    [data swifParams]= readSwif(fullfile(acqPath,swifFiles(end).name),fullfile(acqPath,dimFiles(end).name));
    bmodedata = db(abs(complex(data.I(:,:,1),data.Q(:,:,1))));
    bmodedata = (bmodedata-30)/50;
    bmodedata = bmodedata(:,2:3:end,:);
    z2 = (0:size(bmodedata,1)-1)*(par.c/1e3)/(2*par.fs);
    x1 = (0:size(dttpsLstl,2)-1)*mean(diff(x1));
    x1 = x1-mean(x1);
    x2 = (0:size(bmodedata,2)-1)*mean(diff(x1));
    x2 = x2-mean(x2);
    x1 = x1(:);
    x2 = x2(:);
    save(dttps_file,'bmodedata','z2','x1','x2','-append');
end
%%
clear dttpsLmtl dttpsRmtl dttpsLstl dttpsRstl

clear data
[swei_mtl swei_stl] = deal(struct('cData',[],'x',[],'z',[],'regressionFilter',regressionFilters,'medianFilter',medianFilters));
swei_mtl.x = single(x);
swei_mtl.z = single(z1);
swei_stl.x = single(x);
swei_stl.z = single(z1);
arfi_inv.x = single(x);
arfi_inv.z = single(z);
arfi_inv.regressionFilter = 5+[0:length(regressionFilters)-1]; % TIMESTEPS not FILTER SIZES
b_mode.z = single(z2);
b_mode.x = single(x2);
b_mode.regressionFilter = 0;

%%
N = length(regressionFilters)*length(medianFilters);
ii = 0;
if usejava('jvm');
    H = waitbar(0,'Smoothing...');
else
    fprintf('Smoothing...');
end

ttps(repmat(abs(dlat)>maxDX,[size(ttps,1) 1 1])) = nan;
ttps1(repmat(abs(dlat)>maxDX,[size(ttps1,1) 1 1])) = nan;

[Z X X1] = ndgrid(z1,x,x);
MSK = (X-X1)>0;
clear Z X X1
ttpsL = ttps1;
ttpsR = ttps1;
ttpsR(MSK) = nan;
ttpsL(~MSK) = nan;
offset = round(min(length(x),max(0.5,0.5*(abs(z1-22)))));
offset(:)=0;  %no offset until i figure how it is set;
ttpsL1 = ttpsL;
ttpsR1 = ttpsR;
fprintf('shifting STL...')
for i = 1:length(z1);
    ttpsL1(i,:,:) = interp1(1:length(x),squeeze(ttpsL(i,:,:)).',(1:length(x))-offset(i),'linear').';
    ttpsR1(i,:,:) = interp1(1:length(x),squeeze(ttpsR(i,:,:)).',(1:length(x))+offset(i),'linear').';
end
fprintf('done\n');
ttpsL1(~MSK) = ttpsR1(~MSK);
ttps1 = ttpsL1;
clear ttpsL* ttpsR*


%         dttpsMTL = -1*squeeze(nanmedian(diff(ttps,1,2),3));
%         dttpsSTL = squeeze(nanmedian(diff(ttps1,1,3),2));
%         ttpsL = ttps;
%         ttpsR = ttps;
%         ttpsL(repmat(dlat>0,[size(ttps,1) 1 1])) = nan;
%         ttpsR(repmat(dlat<0,[size(ttps,1) 1 1])) = nan;
%         dttpsMTLR = -1*squeeze(nanmedian(diff(ttpsR,1,2),3));
%         dttpsSTLR = squeeze(nanmedian(diff(ttpsR,1,3),2));
%         dttpsMTLL = -1*squeeze(nanmedian(diff(ttpsL,1,2),3));
%         dttpsSTLL = squeeze(nanmedian(diff(ttpsL,1,3),2));



for fidx = 1:length(regressionFilters);
    filtSize = regressionFilters(fidx);
    %% Loop over regression kernels
    %             if filtSize == 0
    %                 dttpsMTLfilt = dttpsMTL;
    %                 dttpsSTLfilt = dttpsSTL;
    %             else
    %                 dttpsMTLfilt = smoothdim(dttpsMTL,filtSize,'lowess',2);
    %                 dttpsSTLfilt = smoothdim(dttpsSTL,filtSize,'lowess',2);
    %             end
%     MTL=dx./linreg(-1*ttps,max(2,filtSize),2);
%     MTL(MTL>10)=nan;
%     STL=dx./linreg(permute(ttps1,[1 3 2]),max(2,filtSize),2);
%     STL(STL>10)=nan;
%     dxdtMTL = nanmedian(MTL,3);
%     dxdtSTL = nanmedian(STL,3);
            dxdtMTL = dx./nanmedian(linreg(-1*ttps,max(2,filtSize),2),3);
            dxdtSTL = dx./nanmedian(linreg(permute(ttps1,[1 3 2]),max(2,filtSize),2),3);

    
    if D2filt==1
        for modeidx = 1:length(modalities)
            switch(modalities{modeidx})
                case 'SWEI_MTL';
                    %swei_mtl.cData(:,:,fidx,midx) = single(medfilt2((dx./dttpsMTLfilt).^2,[mfiltSize mfiltSize]));
                    swei_mtl.cData(:,:,1,1) = single(medfilt2(dxdtMTL.^2,[medianFilters(1) medianFilters(2)]));
                case 'SWEI_STL'
                    %swei_stl.cData(:,:,fidx,midx) = single(medfilt2((dx./dttpsSTLfilt).^2,[mfiltSize mfiltSize]));
                    swei_stl.cData(:,:,1,1) = single(medfilt2((dxdtSTL).^2,[medianFilters(1) medianFilters(2)]));
            end
        end
    else
        
        for midx = 1:length(medianFilters)
            %% Loop over median filters
            ii = ii+1;
            if usejava('jvm')
                waitbar(ii/N,H);
            end
            mfiltSize = medianFilters(midx);
            for modeidx = 1:length(modalities)
                switch(modalities{modeidx})
                    case 'SWEI_MTL';
                        %swei_mtl.cData(:,:,fidx,midx) = single(medfilt2((dx./dttpsMTLfilt).^2,[mfiltSize mfiltSize]));
                        swei_mtl.cData(:,:,fidx,midx) = single(medfilt2(dxdtMTL.^2,[mfiltSize mfiltSize]));
                    case 'SWEI_STL'
                        %swei_stl.cData(:,:,fidx,midx) = single(medfilt2((dx./dttpsSTLfilt).^2,[mfiltSize mfiltSize]));
                        swei_stl.cData(:,:,fidx,midx) = single(medfilt2((dxdtSTL).^2,[mfiltSize mfiltSize]));
                end
            end
        end
    end
end
if usejava('jvm')
    close(H);
else
    fprintf('done\n');
end
swei_mtl.cData = single(swei_mtl.cData);
swei_stl.cData = single(swei_stl.cData);
%        arfi_inv.cData = single(arfi_inv.cData);
%        b_mode.cData = single(b_mode.cData);

imData = struct('swei_mtl',swei_mtl,'swei_stl',swei_stl);
