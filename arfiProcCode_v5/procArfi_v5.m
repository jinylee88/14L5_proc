function Resfiles = procArfi_v5(varargin)
tic

options.debug = struct(...
    'keyboardFrame',0 ...
    ,'keyboardEnd',0 ...
    ,'overwrite',0 ...
    );

options.save = struct(...
    'Frames',1 ...
    ,'All',0 ...
    ,'compress',1 ...
    );

options.display = struct(...
    'enable',0 ...
    ,'arfiLimit',[0 10] ...
    ,'dxdtLimit',[0 6] ...
    ,'axisLimits',[-10 10 0 30] ...
    ,'arfiFrameTms',0.5 ...
    ,'bmodeMask',0 ...
    ,'depthCorrect',0 ...
    );

options.estimator = struct(...
    'FrameIndices','all' ...
    ,'estimator','PesaventoFlux' ...
    ,'interpFactor',2 ...
    ,'kernelLength',5 ...
    ,'searchRegion',1 ...
    );


options.motionfilter = struct(...
    'enable',1 ...
    ,'order',1 ...
    ,'timeRange',[-inf -0.15] ...
    ,'kickBackFilter',0);

options.bmode = struct(...
    'dynamicRange',35 ...
    ,'offset',10 ...
    ,'scanconvert', 0 ...
    );

options.sws = struct(...
    'enable',1 ...
    ,'dxdtRange',[0.5 8] ...
    ,'xRange',[0.25 5] ...
    ,'tRange',[0 inf] ...
    ,'kernel',[5 1] ...
    );
options.sws.method = {'DTTPS','SMURF'};

options.display.figureHandles = struct(...
    'BMode',100,...
    'IQBMode',101,...
    'ARFI',102,...
    'SWEI',[103:106],...
    'rtSlices',107,...
    'rtMov',108,...
    'SMURF',109,...
    'MMODE',110 ...
    );

if mod(nargin,2)==0
    fNames = dir(fullfile(pwd,'SWIF_AData*.bin'));
    filenames = fullfile(pwd,'data',fNames(end).name);
    args = varargin;
    clear fNames
else
    filenames = varargin{1};
    if nargin>1
        args = varargin(2:end);
    else
        args = {};
    end
end


optfields = fieldnames(options);
for i = 1:2:length(args)
    for idx = 1:length(optfields)
        Tag = optfields{idx};
        if strcmpi(args{i}(1:min(length(Tag),length(args{i}))),Tag)
            subarg = args{i}(length(Tag)+1:end);
            Tagopts = fieldnames(options.(Tag));
            if ~any(strcmpi(Tagopts,subarg))
                warning(sprintf('Could not find option %s.%s. Skipping this option.',Tag,subarg))
            else
                newarg = Tagopts{strcmpi(Tagopts,subarg)};
                options.(Tag).(newarg) = args{i+1};
            end
        end
    end
end
clear Tag Tagopts args idx newarg optfields subarg i

if iscell(filenames)
    fileindices = 1:length(filenames);
else
    fileindices = 1:size(filenames,1);
end

if ispc
    SCPath = 'C:\pjh7\SC2000';
else
    SCPath = '/getlab/pjh7/SC2000';
end
toc1 = toc;
fprintf('Setting up Paths...')
s = path;
if ~any(strfind(s,fullfile(SCPath,'arfiProcCode_v5')));
    addpath(fullfile(SCPath,'arfiProcCode_v5'));
end
if ~any(strfind(s,fullfile(SCPath,'scanconversion')));
    addpath(fullfile(SCPath,'scanconversion'));
end
set(0,'defaultaxesfontweight','bold')

tmpdir = fullfile(SCPath,'tmp');
fprintf(' done (%0.2fs)\n',toc-toc1);
clear s

resfile = '';
Resfiles = [];
for fileidx = fileindices
    

    
    if iscell(filenames)
        fname = filenames{fileidx};
    else
        fname = filenames(fileidx,:);
    end
    [datadir] = fileparts(fname);
    if isempty(datadir)
        datadir = pwd;
    end
    
    % Pull out time stamp from filename
    timeStamp = RetrieveTimeStamp(fname);
    
        if options.save.All
        [pth] = fileparts(fname);
        outputdatadir = fullfile(pth,'res',options.estimator.estimator);
        if ~exist(outputdatadir,'dir')
            mkdir(outputdatadir)
        end
        resfile_all = fullfile(outputdatadir,sprintf('res_all_%s.mat',timeStamp));
        
    if ~options.debug.overwrite && exist(resfile_all,'file')
        fprintf('%s found. Skipping.',resfile_all)
        resfile = resfile_all;
        continue
    end
            
        
        end
    
    % Check to make sure data, dimensions, and parameters files all exist
    dimsname = fullfile(datadir,['SWIF_ADataDims_' timeStamp '.txt']);
    if ~exist(fname, 'file'),error('No data file found');end
    if ~exist(dimsname, 'file'),error('No dimensions file found');end
    
    toc1 = toc;
    % Pull out IQ data
    fprintf('reading SWIF data...')
    [data swifParams]= readSwif(fname, dimsname);
    I = single(data.I);
    Q = single(data.Q);
    
    if strcmpi(options.estimator.interpFactor,'auto')
        options.estimator.interpFactor = floor(1000/size(I,1));
    end
    
    clear data
    fprintf('done (%0.2fs)\n',toc-toc1);
    
    fprintf('Loading parameters...');toc1 = toc;
    parFile = fullfile(datadir,sprintf('par_%s.mat',timeStamp));
    if exist(parFile,'file')
        par = load(parFile);
    end
    
    if ~isfield(swifParams,'imagingMode')
        switch  swifParams.LinesPerSlice
            case 240
                swifParams.imagingMode = 'arfiswei';
            case {120, 256}
                swifParams.imagingMode = 'arfi';
            case 242
                swifParams.imagingMode = 'bmode';
            case {640,320,352,208,416}
                swifParams.imagingMode = 'mmode';
            case 16
                swifParams.imagingMode = 'swei';
            otherwise
                %swifParams.imagingMode = 'mmode';
                %swifParams.imagingMode = input('Input Imaging Mode:','s');
                error('Unable to Parse Mode');
        end
    end
    if strcmpi(swifParams.imagingMode,'arfiL')
        swifParams.imagingMode = 'arfiswei_L';
    end
    
    if any(strfind(swifParams.imagingMode,'_'));
        swifParams.imagingModeShort = swifParams.imagingMode(1:strfind(swifParams.imagingMode,'_')-1);
    else
        swifParams.imagingModeShort = swifParams.imagingMode;
    end
    
    fprintf('Imaging mode is %s\n',swifParams.imagingMode)
    
    if strcmp(swifParams.imagingMode,'bmode')
        showBmode

    else
        
        
        
        parFile1 = fullfile(tmpdir,swifParams.imagingModeShort,'parameters.mat');
        if ~exist(parFile,'file')
            fprintf('Warning: %s not found. Using %s instead.\n',parFile,parFile1)
            parFile = parFile1;
            par = load(parFile);
        end
        
        % Load parameters file
        if exist(parFile1,'file')
            par1 = load(parFile1);
            
            %%%%%%%%%%%%%% TEMPORARY FIX %%%%%%%%%%%%%%%
            if ~strcmpi(swifParams.imagingModeShort,'arfi')
                par.trackParams.paramList.txXyzGridParams = par1.trackParams.paramList.txXyzGridParams;
                par.trackParams.paramList.rxXyzDelayGridParams = par1.trackParams.paramList.rxXyzDelayGridParams;
            else
                par1.trackParams.paramList.txXyzGridParams = par.trackParams.paramList.txXyzGridParams;
                par1.trackParams.paramList.rxXyzDelayGridParams = par.trackParams.paramList.rxXyzDelayGridParams;
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        fprintf('done (%0.2fs)\n',toc-toc1);
        if ~isfield(par,'pushbeamNum')
            par.pushbeamNum = ones(1,par.numBeamGroups)*((par.nBeams+1)/2);
        end
        
        clear SCPath datadir dimsname parFile parFile2
        
        FsMhz = par.fs; % sampling frequency per beam
        nParRx = par.nBeams; % number of parallel receive beams
        par.interpFactor = options.estimator.interpFactor;
        if ~isfield(par,'beamNum')
            par.beamNum = 1:par.numBeamGroups;
        end
        % Compute axial vector
        c = 1.540/2; % speed of sound (mm/us)
        N = size(I,1)*options.estimator.interpFactor;
        axial0 = (0:N-1)*c/(FsMhz*options.estimator.interpFactor);
        
        % generate time vector
        % time = 0 is the last push for the multi-focal zone sequences
        if length(par.priusec) == 1
            if ~isfield(swifParams,'numTs')
                swifParams.numTs = par.ensemble;
            end
            t = ((-par.npush*length(par.push_focaldepth)-par.nref+1):(swifParams.numTs-par.nref-par.npush*length(par.push_focaldepth)))*par.priusec*1e-3;
        else
            priUsec = [max(par.priusec)*ones(1,par.nref+length(par.pushFocalDepth)*par.npush)];
            for i = 1:length(par.priusec);
                priUsec = [priUsec par.priusec(i)*ones(1,par.ntrack(i))];
            end
            priUsec = [-1*sum(priUsec(1:par.nref+length(par.push_focaldepth)*par.npush-1)) priUsec(1:end-1)];
            t = cumsum(1e-3 * priUsec);
        end
        
        %Calculate frame count
        if ~isfield(swifParams,'numFrames')
            par.numFrames = max(1,floor(swifParams.FrameCount./length(t)));
        else
            par.numFrames = swifParams.numFrames;
        end
        if mod(par.numFrames,1)
            if strcmpi(swifParams.imagingModeShort,'swei') %adjust if necessary
                t = t(1:50);
                par.numFrames = swifParams.FrameCount./length(t);
            end
        end
        
        %sws.dxdt = zeros(length(sws.z),par.numFrames);
        %sws.resid = sws.dxdt;
        
        if strcmpi(swifParams.imagingModeShort,'mmode');
            
            I = reshape(permute(reshape(I,size(I,1),par.nBeams,par.numBeamGroups,par.ensemble,par.numAcq),[1 2 4 3 5]),size(I,1),par.nBeams,[]);
            Q = reshape(permute(reshape(Q,size(Q,1),par.nBeams,par.numBeamGroups,par.ensemble,par.numAcq),[1 2 4 3 5]),size(Q,1),par.nBeams,[]);
            par.numFrames = par.numFrames*par.numBeamGroups;
        end
        
        if ~isfield(options.estimator,'FrameIndices') || strcmpi(options.estimator.FrameIndices,'all')
            Frameindices = 1:par.numFrames;
        else
            Frameindices = options.estimator.FrameIndices(options.estimator.FrameIndices<=par.numFrames);
        end
        clear bmodedata_all bimg_all arfidata_all cdata sws_all cc_all arfidata0_all
        for frmii = 1:length(Frameindices)
            frmidx = Frameindices(frmii);
            %Compute Displacements
            toc1 = toc;
            fprintf('Computing Frame %d of %d...',frmidx,par.numFrames);
            idx = (1:length(t))+(frmidx-1)*length(t);
            switch lower(options.estimator.estimator)
                case 'loupas'
                    [arfidata cc_coef Iup Qup] = runLoupas(I(:,:,idx),Q(:,:,idx),options.estimator.interpFactor, options.estimator.kernelLength, axial0, par);
                case 'kasai'
                    [arfidata cc_coef Iup Qup] = runKasai(I(:,:,idx),Q(:,:,idx),options.estimator.interpFactor, options.estimator.kernelLength, axial0, par);
                case 'flux'
                    [arfidata cc_coef Iup Qup] = runFlux(I(:,:,idx),Q(:,:,idx),t,options.estimator.interpFactor, options.estimator.kernelLength, axial0, par);
                case 'loupasflux'
                    [arfidata cc_coef Iup Qup] = runLoupasFlux(I(:,:,idx),Q(:,:,idx),t,options.estimator.interpFactor, options.estimator.kernelLength, axial0, par);
                case 'pesavento'
                    [arfidata cc_coef Iup Qup] = runPesavento(I(:,:,idx),Q(:,:,idx),options.estimator.interpFactor,options.estimator.kernelLength*2,options.estimator.searchRegion,axial0, par);
                case 'pesaventoflux'
                    [arfidata cc_coef Iup Qup] = runPesaventoFlux(I(:,:,idx),Q(:,:,idx),t,options.estimator.interpFactor,options.estimator.kernelLength*2,options.estimator.searchRegion,axial0, par);
                otherwise
                    error('incompatible Displacement Method ''%s''',options.estimator.estimator);
            end
            fprintf('done (%0.2fs)\n',toc-toc1);
            
            if options.motionfilter.kickBackFilter
                toc1 = toc;
                fprintf('Calculating Transducer Kickback...');
                [kickback] = removeKickBack(I(:,:,idx),Q(:,:,idx),options.estimator.interpFactor,options.estimator.kernelLength,options.estimator.searchRegion,axial0, par);
                arfidata = arfidata-kickback(1:size(arfidata,1),:,:);
                fprintf('done (%0.2fs)\n',toc-toc1);
            end
            
            
            arfidata = single(arfidata);
            Iup = single(Iup);
            Qup = single(Qup);
            
            nPushLocations = size(arfidata,2)/nParRx; % Compute the number of push locations actually used in acquisition
            
            % Remove extra axial samples
            axial = axial0(1:size(arfidata,1));
            [lat] = genLatMatrix(par);
            apex = par.trackParams.specList.scanSpecs.apexMm(3);
            theta = lat(1,:);
            
            
            arfidata0= arfidata;
            if options.motionfilter.enable
                toc1 = toc;
                fprintf('Motion filtering...')
                tmask = false(size(t));
                for i = 1:2:length(options.motionfilter.timeRange)
                    tmask(t>options.motionfilter.timeRange(i) & t<options.motionfilter.timeRange(i+1)) = true;
                end
                tmask(par.nref+(1:par.npush*length(par.pushFocalDepth))) = false;
                arfidata = linearmotionfilter(arfidata0,t,find(tmask),options.motionfilter.order);
                fprintf('done (%0.2fs)\n',toc-toc1);
                
            end
            
            %Cubic Interpolate Push and Reverb
            tidx1 = [par.nref+[-1 0] par.nref+par.npush+[2:3]];
            tidx2 = [par.nref+[1:par.npush+1]];
            [residtmp motion1] = linearmotionfilter(arfidata0,t,tidx1,3);
            arfidata0(:,:,tidx2) = motion1(:,:,tidx2);
            [residtmp motion1] = linearmotionfilter(arfidata,t,tidx1,3);
            arfidata(:,:,tidx2) = motion1(:,:,tidx2);
            clear residtmp motion1 tidx1 tidx2;
            
            
            cc_coef = reshape(cc_coef, [size(arfidata,1) nParRx nPushLocations size(arfidata,3)]); % reshape into 4-D matrix
            arfidata = reshape(arfidata, [size(arfidata,1) nParRx nPushLocations size(arfidata,3)]); % reshape into 4-D matrix
            arfidata0 = reshape(arfidata0, [size(arfidata0,1) nParRx nPushLocations size(arfidata0,3)]); % reshape into 4-D matrix
            
            % Trim lateral vector to correct number of push locations if necessary
            if par.numBeamGroups ~= nPushLocations
                if ~strcmp(swifParams.imagingModeShort,'mmode')
                    lat = lat((par.numBeamGroups+1)/2 + [-floor(nPushLocations/2):floor(nPushLocations/2)],:);
                end
            end
            
            toc1 = toc;
            fprintf('Extracting B-Mode...')
            
            par.imagingdepth = 0.1*max(axial0);
            IQ = complex(Iup(1:size(arfidata,1),:,:),Qup(1:size(arfidata,1),:,:));
            IQ = reshape(IQ, [size(arfidata,1) nParRx nPushLocations length(t)]); % reshape into 4-D matrix
            
            if any(max(abs(diff(lat)))>0.1)
                blat = linspace(min(lat(:)),max(lat(:)),min(512,length(unique(round(1e4*lat(:))))));
                blat = unique(round(lat*1e4))/1e4;
                IQsws = [];
                if strcmpi(swifParams.imagingMode,'mmode')
                    IQsws = interp1(lat(mod(frmidx-1,par.numBeamGroups)+1,:),squeeze(mean(IQ(:,:,1,1:par.nref),4))',blat)';
                    cc1  = interp1(lat(mod(frmidx-1,par.numBeamGroups)+1,:),squeeze(1/255*mean(double(cc_coef(:,:,1,(par.nref+4):end)),4))',blat)';
                    
                else
                    for i = 1:size(arfidata,3)
                        IQsws(:,:,i) = interp1(1e-4*round(lat(i,:)*1e4),squeeze(mean(IQ(:,:,i,1:par.nref),4))',blat)';
                        %if flags.cc
                            cc1(:,:,i) = interp1(1e-4*round(lat(i,:)*1e4),squeeze(1/255*mean(double(cc_coef(:,:,i,(par.nref+4):end)),4))',blat)';
                        %end
                    end
                end
                bmodedata0 = nanmean(abs(IQsws),3);
                interpidx = all(all(isnan(IQsws),3));
                bmodedata0(:,interpidx) = interp1(blat(~interpidx),bmodedata0(:,~interpidx)',blat(interpidx))';
            else
                blat = mean(lat,1);
                bmodedata0 = mean(mean(abs(IQ(:,:,:,1:par.nref)),4),3);
                %bmodedata0 = reshape(permute(abs(IQ),[1 2 4 3]),size(IQ,1),size(IQ,2),[]);
            end
            iq0 = min(1,max(0,(db(IQ)+options.bmode.offset-options.bmode.dynamicRange)/options.bmode.dynamicRange));
            bmodedata0 = min(1,max(0,(db(bmodedata0)+options.bmode.offset-options.bmode.dynamicRange)/options.bmode.dynamicRange));
            if ~exist('bmodedata_all','var')
                bmodedata_all = uint8(zeros([size(bmodedata0) length(Frameindices)]));
            end
            bmodedata_all(:,:,frmii) = uint8(255*bmodedata0);
            
            scB = struct(...
                'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
                'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
                'latinc',.1e-3,...
                'axialmin',-1e-3*apex,...
                'axialmax',1e-2*(par.imagingdepth) - 1e-3*apex,...
                'axialinc',.1e-3,...
                'min_phi',min(blat),...
                'span_phi',max(blat)-min(blat),...
                'apex',0.1*apex,...
                'fsiq',par.fs*options.estimator.interpFactor*1e6...
                );
            if options.bmode.scanconvert
                [bimg baxial blat] = Scan_Convert_Sector(single(bmodedata0),scB);
                bimg_all(:,:,frmii) = bimg;
            end
            
            fprintf(' done (%0.2fs)\n',toc-toc1);
            
            
            
            switch lower(swifParams.imagingModeShort)
                case 'swei'
                    %arfidata = squeeze(arfidata);
                    %cc = squeeze(cc);
                    bmwidth = 2 * 1e3 * (par.c/(par.F0MHz*1e6))' * par.pushFnum;
                    xfoc = sind(lat(mod(frmidx-1,par.numBeamGroups)+1,:))*par.pushFocalDepth(1);
                    Ridx = find(xfoc>0*bmwidth);
                    Lidx = find(xfoc<0*(-bmwidth));
                    Cidx = [];%find(xfoc>(-bmwidth)&xfoc<bmwidth);
                    Cidx1 = 1:size(lat,2);
                case 'mmode'
                    sweidata0 = permute(arfidata0,[1 2 4 3]);
                    sweidata = permute(arfidata,[1 2 4 3]);
                    [lat dlat] = genLatMatrix(par);
                    theta = lat(1,:);
                    dtheta = dlat(1,:);
                    bmwidth = 2 * 1e3 * (par.c/(par.F0MHz*1e6))' * par.pushFnum;
                    xfoc = sind(lat(1,:))*par.pushFocalDepth(1);
                    Ridx = find(xfoc>0*bmwidth);
                    Lidx = find(xfoc<0*(-bmwidth));
                    Cidx = [];%find(xfoc>(-bmwidth)&xfoc<bmwidth);
                    Cidx1 = 1:size(lat,2);
                case 'arfi'
                    if ~strcmpi(par.mode,'mmode')
                        arfidata = reshape(arfidata,length(axial),par.numBeamGroups*par.nBeams,length(t));
                        arfidata0 = reshape(arfidata0,length(axial),par.numBeamGroups*par.nBeams,length(t));
                        cc_coef = reshape(cc_coef,length(axial),par.numBeamGroups*par.nBeams,length(t));
                        cc_onax = cc_coef;
                        lat = reshape(lat',1,[]);
                        theta = lat;
                    else
                        arfidata1 = arfidata;
                        [arfidata0 sweidata0 theta dtheta] = splitARFISWEI(arfidata0,par);
                        [arfidata sweidata theta dtheta] = splitARFISWEI(arfidata,par);
                        %[arfidata0 sweidata0 theta dtheta] = splitARFISWEI(arfidata0,par);
                        [cc_onax cc_coef] = splitARFISWEI(cc_coef,par);
                        cc_onax(isnan(cc_onax)) = 0;
                        cc_coef(isnan(cc_coef)) = 0;
                        W = exp(-100*(1-abs(cc_onax)));
                        K = repmat([0.1 0.8 0.1],5,1);
                        arfidata2 = arfidata;
                        arfidata2(isnan(arfidata2)) = 0;
                        arfidata = convn(arfidata2.*W,K,'same')./convn(W.*~isnan(arfidata),K,'same');
                        sweidata2 = sweidata;
                        sweidata2(isnan(sweidata2)) = 0;
                        W = exp(-100*(1-abs(cc_coef)));
                        sweidata = convn(sweidata2.*W,K,'same')./convn(W.*~isnan(sweidata),K,'same');
                        
                    end
                    Cidx = 1:par.nBeams;
                    Cidx1 = Cidx;
                case 'arfiswei'
                    [cc_onax cc_coef] = splitARFISWEI_Centered(cc_coef,par);
                    [arfidata0 sweidata0 theta dtheta] = splitARFISWEI_Centered(arfidata0,par);
                    [arfidata sweidata theta dtheta] = splitARFISWEI_Centered(arfidata,par);
                    bmwidth = 2 * 1e3 * (par.c/(par.F0MHz*1e6))' * par.pushFnum;
                    xfoc = sind(lat(round((par.numBeamGroups+1)/2),:))*max(par.pushFocalDepth);
                    Ridx = 11:16;%find(xfoc>bmwidth);
                    Lidx = 1:6;%find(xfoc<(-bmwidth));
                    Cidx = 7:10;%find(xfoc>(-bmwidth)&xfoc<bmwidth);
                    Cidx1 = 1:size(lat,2);
            end
            
            [TH R] = meshgrid(double(theta),double(axial));
            X = R.*sind(TH) - apex*tand(TH);
            Z = R.*cosd(TH);
            
            
            
            if options.sws.enable && ~strcmpi(par.mode,'arfi');
                toc1 = toc;
                fprintf('Calculating Shear Wave Velocities...')
                [TTP DX vTest] = TTPS(theta,dtheta,axial,apex,t,sweidata,options.sws);
                [ii jj] = sort(par.pushbeamNum);
                TTPsort = stdFilt(TTP(:,:,jj),[5 3],1);
                DXsort = DX(:,:,jj);
                if ~iscell(options.sws.method);
                    swsmethod = {options.sws.method};
                else
                    swsmethod = options.sws.method;
                end
                for i = 1:length(swsmethod)
                    fprintf('%s...',options.sws.method{i});
                    switch upper(options.sws.method{i})
                        case 'TTPS'
                            CT_TTPS = sws_TTPS(TTPsort,DXsort);
                        case 'DTTPS'
                            CT_DTTPS = sws_DTTPS(TTPsort,DXsort,1);
                        case 'SMURF'
                            CT_SMURF = permute(sws_DTTPS(permute(TTPsort,[1 3 2]),permute(DXsort,[1 3 2]),1),[1 3 2]);
                    end
                end
                eval(sprintf('CT = CT_%s;',upper(swsmethod{1})));
                fprintf('done (%0.2fs)\n',toc-toc1);
                sws_all(:,:,frmii) = nanmedian(CT,3);
            end
            
            sc = struct(...
                'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
                'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
                'latinc',.1e-3,...
                'axialmin',-1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
                'axialmax',1e-2*(par.imagingdepth)- 1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
                'axialinc',.1e-3,...
                'min_phi',min(lat(:)),...
                'span_phi',max(lat(:))-min(lat(:)),...
                'apex',0.1*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
                'fsiq',par.fs*options.estimator.interpFactor*1e6...
                );
            
            if options.display.enable
                toc1 = toc;
                fprintf('Generating Display...')
                realtime_display_v5
                fprintf('done (%0.2fs)\n',toc-toc1);
            end
            
            
            
            
            if options.save.Frames
                %outputdatadir = fullfile('C:\pjh7\SC2000\RES\',datestr(now,'yyyymmdd'));
                [pth] = fileparts(fname);
                outputdatadir = fullfile(pth,'res');
                
                % Save time stamped results file
                if par.numFrames>1
                    outputdatadir = fullfile(pth,'res',options.estimator.estimator,sprintf('res_%s',timeStamp));
                    resfile(frmidx,:) = fullfile(outputdatadir,sprintf('res_%s_%03.0f.mat',timeStamp,frmidx));
                else
                    outputdatadir = fullfile(pth,'res',options.estimator.estimator);
                    resfile = fullfile(outputdatadir,sprintf('res_%s.mat',timeStamp));
                end
                
                
                if ~exist(outputdatadir,'dir');mkdir(outputdatadir);end
                fprintf('Saving %s...\n',resfile(frmidx,:))
                cc = uint8(min(255,512*max(0,cc_coef-0.5)));
                
                if options.save.compress
                arfi_scale = 255/2^15;
                ccmap = 1-10.^(-6*(0:255)/255);
                arfidata = int16(arfidata/arfi_scale);
                arfidata0 = int16(arfidata0/arfi_scale);
                cc = uint8(log10(1-abs(squeeze(cc_coef)))*(-255/6));
                bmodedata0 = uint8(bmodedata0*255);
                iq0 = uint8(iq0*255);
                if exist('sweidata','var')
                    sweidata = int16(sweidata/arfi_scale);
                    sweidata0 = int16(sweidata0/arfi_scale);
                end
                if exist('vTest','var')
                    vTest = int16(vTest/arfi_scale);
                end
                switch par.probe        
                    case 'linear'
                        save(resfile(frmidx,:),'arfidata','sweidata','sweidata0','cc','axial','lat','t','bmodedata0','blat','arfidata0','arfi_scale','ccmap','iq0')
                    case 'phased'
                        if exist('sweidata','var')
                            if exist('vTest','var')
                            save(resfile(frmidx,:),'arfidata','sweidata','vTest','cc','axial','theta','apex','t','bmodedata0','blat','sc','scB','TTP','CT_*','arfidata0','sweidata0','arfi_scale','ccmap','lat')
                            else
                            save(resfile(frmidx,:),'arfidata','sweidata','cc','axial','theta','apex','t','bmodedata0','blat','sc','scB','TTP','CT_*','arfidata0','sweidata0','arfi_scale','ccmap','lat')
                            end
                        else
                            save(resfile(frmidx,:), 'arfidata','cc','axial','theta','apex','t','bmodedata0','blat','sc','scB','arfidata0','arfi_scale','ccmap','lat');
                        end
                end
                
                else
                if exist('sweidata','var')
                    save(resfile(frmidx,:),'arfidata','sweidata','cc','axial','theta','apex','t','bmodedata0','btheta','sc','scB','TTP','CT_*','arfidata0')
                    
                else
                    save(resfile(frmidx,:), 'arfidata','cc','axial','theta','apex','t','bmodedata0','btheta','sc','scB','arfidata0');
                end
                end
                toc
            end
            
            if ~exist('arfidata_all','var')
                arfidata_all = int16(zeros([size(squeeze(arfidata)) length(Frameindices)]));
                arfidata0_all = arfidata_all;
                cc_all = uint8(arfidata_all);
            end
            
            
            if length(size(squeeze(arfidata))) == 4
            arfidata_all(:,:,:,:,frmii) = int16(squeeze(arfidata)/arfi_scale);
            arfidata0_all(:,:,:,:,frmii) = int16(squeeze(arfidata0)/arfi_scale);
            cc_all(:,:,:,:,frmii) = uint8(log10(1-abs(squeeze(cc_coef)))*(-255/6));
            else
            arfidata_all(:,:,:,frmii) = int16(squeeze(arfidata)/arfi_scale);
            arfidata0_all(:,:,:,frmii) = int16(squeeze(arfidata0)/arfi_scale);
            cc_all(:,:,:,frmii) = uint8(log10(1-abs(squeeze(cc_onax)))*(-255/6));
            end
%             ARFIDATA0(:,:,:,frmii) = squeeze(arfidata0);
%             
%             if exist('kickback','var')
%             KICKBACK(:,:,:,frmii) = squeeze(kickback);
%             end
%             
            
            if options.debug.keyboardFrame
                keyboard
            end
        end
        ii = reshape(reshape(1:32,[],2)',[],1);

        arfi_scale = 255/2^15;
        ccmap = 1-10.^(-6*(0:255)/255);

        
        if options.save.All
            fprintf('saving %s...',resfile_all)
            save(resfile_all, 'arfidata_all','arfidata0_all','cc_all','axial','theta','apex','t','bmodedata_all','btheta','sc','scB','ccmap','arfi_scale');
            fprintf('done\n');
            resfile = resfile_all
        end
        
    end
    
    if isempty(Resfiles)
        Resfiles = resfile;
    else
        Resfiles = char(Resfiles,resfile);
    end
    
    
    if options.debug.keyboardEnd
        keyboard
    end
end
end


function [arfidata sweidata lat dlat] = splitARFISWEI(arfidata0,par)
arfidata = single(zeros(size(arfidata0,1),par.numBeamGroups,par.ensemble));
sweidata = permute(reshape(arfidata0,size(arfidata0,1),par.nBeams,par.numBeamGroups,par.ensemble),[1 2 4 3]);
[lat dlat] = genLatMatrix(par);
dii = [-1 0 1];
w = [0.25 0.5 0.25];
W = repmat(w,[size(sweidata,1),1,size(sweidata,3)]);
for i = 1:par.numBeamGroups;
    idx0 = max(1,min(size(sweidata,2),floor(par.pushbeamNum(i))+dii));
    idx1 = max(1,min(size(sweidata,2),ceil(par.pushbeamNum(i))+dii));
    frac = mod(par.pushbeamNum(i),1);
    arfidata(:,par.pushbeamNum(i),:) = ...
        (1-frac)*sum(sweidata(:,idx0,:,i).*W,2)./sum(W,2) +...
        (frac)*sum(sweidata(:,idx1,:,i).*W,2)./sum(W,2);
end
lat = lat(1,:);
end

function [arfidata sweidata lat dlat] = splitARFISWEI_Centered(arfidata0,par)
arfidata = single(zeros(size(arfidata0,1),par.numBeamGroups,par.ensemble));
sweidata = permute(reshape(arfidata0,size(arfidata0,1),par.nBeams,par.numBeamGroups,par.ensemble),[1 2 4 3]);
[lat dlat] = genLatMatrix(par);
dii = [-1 0 1];
w = [0.25 0.5 0.25];
W = repmat(w,[size(sweidata,1),1,size(sweidata,3)]);
for i = par.beamNum;
    idx0 = max(1,min(size(sweidata,2),floor(par.pushbeamNum(i))+dii));
    idx1 = max(1,min(size(sweidata,2),ceil(par.pushbeamNum(i))+dii));
    frac = mod(par.pushbeamNum(i),1);
    arfidata(:,i,:) = ...
        (1-frac)*nansum(sweidata(:,idx0,:,i).*W,2)./nansum(~isnan(sweidata(:,idx0,:,i)).*W,2) +...
        (frac)*nansum(sweidata(:,idx1,:,i).*W,2)./nansum(~isnan(sweidata(:,idx1,:,i)).*W,2);
end
lat = lat(1,:);
end





function [lat dlat] = genLatMatrix(par)
txBeams = interp1(1:par.numBeamGroups,par.trackParams.paramList.txXyzGridParams.P,par.beamNum(:));
rxBeams = par.trackParams.paramList.rxXyzDelayGridParams.beamDeltaP;
lat = repmat(txBeams(:),[1 length(rxBeams)])'+repmat(rxBeams, [length(txBeams) 1])';
switch par.probe
    case 'linear'
        txSpacing = nanmean(squeeze(par.trackParams.paramList.txXyzGridParams.fociMm(1,:,:))./squeeze(par.trackParams.paramList.txXyzGridParams.P)); % nanmean in case there is a beam located at 0
        lat = lat'*txSpacing;
        if isfield(par,'pushbeamNum')
            txPush = interp1(1:par.numBeamGroups,par.trackParams.paramList.txXyzGridParams.P,par.pushbeamNum);
            dlat = lat - repmat(txSpacing*txPush(:),[1 size(lat,2)]);
        else
            dlat = lat;
        end
    case 'phased'
        txSpacing = (par.trackParams.paramList.txXyzGridParams.thetaX_DEG./par.trackParams.paramList.txXyzGridParams.P);
        lat = (lat.*repmat(txSpacing',[size(lat,1) 1]))';
        if isfield(par,'pushbeamNum')
            txPush = interp1(1:par.numBeamGroups,par.trackParams.paramList.txXyzGridParams.thetaX_DEG,par.pushbeamNum);
            dlat = lat - repmat(txPush(:),[1 size(lat,2)]);
        else
            dlat = lat;
        end
end
end
