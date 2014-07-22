function Resfiles = procArfi_Clinical(varargin)
tic

flags = struct(...
    'calcSWS',          1 ...
    ,'saveRes',         0 ...
    ,'displayflag',     1 ...
    ,'keyboard',        1 ...
    );

mf = struct(...
    'Enable',       1,...
    'Order',        1,...
    'tMotion_ms',   [-inf -0.1 ])%2.25 inf]);

rt = struct(...
    'dispLim',     [0 10],...
    'axisLimits',    [-10 10 0 30],...
    't_arfi_ms', 0.5, ...  %0.6
    'bdB',      30,...
    'bdBoffset',0,...
    'bMask',1,...
    'FrameIndices','all',...
    'depthCorrect',0);

dispEst = struct(...
    'method','flux',...
    'interpFactor',1,...
    'kernelLength',5,...
    'searchRegion',1,...
    'removeKickBack',0);

sws = struct(...
    'dxdtrange',    [0.5 8]...
    ,'xrange',       [0.25 5]...
    ,'trange',       [0 inf]...
    ,'smoothingkernel',[5 1]);
sws.method = {'SMURF','DTTPS','TTPS'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figHandles = struct(...
    'BMode',100,...
    'IQBMode',101,...
    'ARFI',102,...
    'SWEI',[103:106],...
    'rtSlices',107,...
    'rtMov',108,...
    'SMURF',109,...
    'MMODE',110 ...
    );

flags.cc = 1;
flag.savesws = 1;
rt.ccRng = [0.95 1];


if mod(nargin,2)==0
    fNames = dir(fullfile(pwd,'SWIF_AData*.bin'));
    filenames = fullfile(pwd,'data',fNames(end-1).name);
    args = varargin;
else
    filenames = varargin{1};
    if nargin>1
    args = varargin(2:end);
    else
    args = {};
    end
end

for i = 1:2:length(args)
    dispEst.(args{i}) = args{i+1};
end

if iscell(filenames)
    fileindices = 1:length(filenames);
else
    fileindices = 1:size(filenames,1);
end

Resfiles = [];
for fileidx = fileindices
    if iscell(filenames)
        fname = filenames{fileidx};
    else
        fname = filenames(fileidx,:);
    end
    
    resfile = '';
    
if ispc
        SCPath = 'C:\pjh7\SC2000';
    else
        SCPath = '/getlab/pjh7/SC2000';
    end
        toc1 = toc;
        fprintf('Setting up Paths...')
        s = path;
        if ~any(strfind(s,fullfile(SCPath,'arfiProcCode_v4')));
            addpath(fullfile(SCPath,'arfiProcCode_v4'));
        end
        if ~any(strfind(s,fullfile(SCPath,'scanconversion')));
            addpath(fullfile(SCPath,'scanconversion'));
        end
        set(0,'defaultaxesfontweight','bold')
        [datadir] = fileparts(fname);
        if isempty(datadir)
            datadir = pwd;
        end
        tmpdir = fullfile(SCPath,'tmp');
        fprintf(' done (%0.2fs)\n',toc-toc1);
    
    % Pull out time stamp from filename
    timeStamp = RetrieveTimeStamp(fname);
    
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
    
    if strcmpi(dispEst.interpFactor,'auto')
        dispEst.interpFactor = floor(1000/size(I,1));
    end
    
    clear data
    fprintf('done (%0.2fs)\n',toc-toc1);
    
    
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
                swifParams.imagingMode = input('Input Imaging Mode:','s');
                %error('Unable to Parse Mode');
        end
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
        
        
        parFile = fullfile(datadir,sprintf('par_%s.mat',timeStamp));
        parFile1 = fullfile(tmpdir,swifParams.imagingModeShort,'parameters.mat');
        if ~exist(parFile,'file')
            fprintf('Warning: %s not found. Using %s instead.\n',parFile,parFile1)
            parFile = parFile1;
        end
        
        fprintf('Loading parameters...');toc1 = toc;
        % Load parameters file
        par = load(parFile);
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
        
        FsMhz = par.fs; % sampling frequency per beam
        nParRx = par.nBeams; % number of parallel receive beams
        par.interpFactor = dispEst.interpFactor;
        if ~isfield(par,'beamNum')
            par.beamNum = 1:par.numBeamGroups;
        end
        % Compute axial vector
        c = 1.540/2; % speed of sound (mm/us)
        N = size(I,1)*dispEst.interpFactor;
        axial0 = (0:N-1)*c/(FsMhz*dispEst.interpFactor);
        
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
        
        if ~isfield(rt,'FrameIndices') || strcmpi(rt.FrameIndices,'all')
            Frameindices = 1:par.numFrames;
        else
            Frameindices = rt.FrameIndices(rt.FrameIndices<=par.numFrames);
        end
        clear Bmodedata0 Bimg ARFIDATA cdata SWS CC
        for frmii = 1:length(Frameindices)
            frmidx = Frameindices(frmii);
            %Compute Displacements
            toc1 = toc;
            fprintf('Computing Frame %d of %d...',frmidx,par.numFrames);
            idx = (1:length(t))+(frmidx-1)*length(t);
            switch lower(dispEst.method)
                case 'loupas'
                    [arfidata cc Iup Qup] = runLoupas(I(:,:,idx),Q(:,:,idx),dispEst.interpFactor, dispEst.kernelLength, axial0, par);
                case 'kasai'
                    [arfidata cc Iup Qup] = runKasai(I(:,:,idx),Q(:,:,idx),dispEst.interpFactor, dispEst.kernelLength, axial0, par);
                case 'flux'
                    [arfidata cc Iup Qup] = runFlux(I(:,:,idx),Q(:,:,idx),t,dispEst.interpFactor, dispEst.kernelLength, axial0, par);
                case 'loupasflux'
                    [arfidata cc Iup Qup] = runLoupasFlux(I(:,:,idx),Q(:,:,idx),t,dispEst.interpFactor, dispEst.kernelLength, axial0, par);
                case 'pesavento'
                    [arfidata cc Iup Qup] = runPesavento(I(:,:,idx),Q(:,:,idx),dispEst.interpFactor,dispEst.kernelLength,dispEst.searchRegion,axial0, par);
                case 'pesaventoflux'
                    [arfidata cc Iup Qup] = runPesaventoFlux(I(:,:,idx),Q(:,:,idx),t,dispEst.interpFactor,dispEst.kernelLength,dispEst.searchRegion,axial0, par);    
                otherwise
                    error('incompatible Displacement Method ''%s''',dispEst.method);
            end
            fprintf('done (%0.2fs)\n',toc-toc1);
            
            if dispEst.removeKickBack
                toc1 = toc;
                fprintf('Calculating Transducer Kickback...');
               [kickback] = removeKickBack(I(:,:,idx),Q(:,:,idx),dispEst.interpFactor,dispEst.kernelLength,dispEst.searchRegion,axial0, par);
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
           
        
            if mf.Enable
                toc1 = toc;
                fprintf('Motion filtering...')
                tmask = false(size(t));
                for i = 1:2:length(mf.tMotion_ms)
                    tmask(t>mf.tMotion_ms(i) & t<mf.tMotion_ms(i+1)) = true;
                end
                arfidata0 = arfidata;
                arfidata = linearmotionfilter(arfidata0,t,find(tmask),mf.Order);
                fprintf('done (%0.2fs)\n',toc-toc1);
            else
                arfidata0 = arfidata;
            end
            
            %Cubic Interpolate Push and Reverb
            tidx1 = [par.nref+[-1 0] par.nref+par.npush+[2:3]];
            tidx2 = [par.nref+[1:par.npush+1]];
            [residtmp motion1 A] = linearmotionfilter(arfidata0,t,tidx1,3);
            arfidata0(:,:,tidx2) = motion1(:,:,tidx2);
            [residtmp motion1 A] = linearmotionfilter(arfidata,t,tidx1,3);
            arfidata(:,:,tidx2) = motion1(:,:,tidx2);
            clear residtmp motion1;
            
            
%             if flags.savesws
%                 zidx = find(axial>0.75*min(par.push_focaldepth) & axial<max(par.push_focaldepth));
%                 slice(:,:,frmidx) = squeeze(nansum(arfidata(zidx,:,:).*max(0,double(cc(zidx,:,:))-0.99).^2)./nansum(max(0,double(cc(zidx,:,:))-0.99).^2));
%             end
%             
            
            cc = reshape(cc, [size(arfidata,1) nParRx nPushLocations size(arfidata,3)]); % reshape into 4-D matrix
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
            
            if any(diff(lat))
                Th1 = linspace(min(lat(:)),max(lat(:)),min(512,length(unique(lat(:)))));
                IQsws = [];
                if strcmpi(swifParams.imagingMode,'mmode')
                    IQsws = interp1(lat(mod(frmidx-1,par.numBeamGroups)+1,:),squeeze(mean(IQ(:,:,1,1:par.nref),4))',Th1)';
                    if flags.cc
                        CC  = interp1(lat(mod(frmidx-1,par.numBeamGroups)+1,:),squeeze(1/255*mean(double(cc(:,:,1,(par.nref+4):end)),4))',Th1)';
                    end
                else
                    for i = 1:size(arfidata,3)
                        IQsws(:,:,i) = interp1(lat(i,:),squeeze(mean(IQ(:,:,i,1:par.nref),4))',Th1)';
                        if flags.cc
                            CC(:,:,i) = interp1(lat(i,:),squeeze(1/255*mean(double(cc(:,:,i,(par.nref+4):end)),4))',Th1)';
                        end
                    end
                end
                bmodedata0 = nanmean(abs(IQsws),3);
                interpidx = all(all(isnan(IQsws),3));
                bmodedata0(:,interpidx) = interp1(Th1(~interpidx),bmodedata0(:,~interpidx)',Th1(interpidx))';
            else
                Th1 = lat(1,:);
                bmodedata0 = mean(abs(IQ(:,:,:,1:par.nref)),4);
                %bmodedata0 = reshape(permute(abs(IQ),[1 2 4 3]),size(IQ,1),size(IQ,2),[]);
            end
            bmodedata0 = min(1,max(0,(db(bmodedata0)+rt.bdBoffset-rt.bdB)/rt.bdB));
            Bmodedata0(:,:,frmii) = bmodedata0;
            scB = struct(...
                'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
                'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
                'latinc',.1e-3,...
                'axialmin',-1e-3*apex,...
                'axialmax',1e-2*(par.imagingdepth) - 1e-3*apex,...
                'axialinc',.1e-3,...
                'min_phi',min(Th1),...
                'span_phi',max(Th1)-min(Th1),...
                'apex',0.1*apex,...
                'fsiq',par.fs*dispEst.interpFactor*1e6...
                );
            [bimg baxial blat] = Scan_Convert_Sector(single(bmodedata0),scB);
            Bimg(:,:,frmii) = bimg;
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
                    cc = reshape(cc,length(axial),par.numBeamGroups*par.nBeams,length(t));
                    lat = reshape(lat',1,[]);
                    theta = lat;
                    else
                    arfidata1 = arfidata;
                    [arfidata sweidata theta dtheta] = splitARFISWEI(arfidata,par);
                    %[arfidata0 sweidata0 theta dtheta] = splitARFISWEI(arfidata0,par);
                    [cc_onax cc] = splitARFISWEI(cc,par);
                    cc_onax(isnan(cc_onax)) = 0;
                    cc(isnan(cc)) = 0;
                    W = exp(-100*(1-abs(cc_onax)));
                    K = repmat([0.1 0.8 0.1],5,1);
                    arfidata2 = arfidata;
                    arfidata2(isnan(arfidata2)) = 0;
                    arfidata = convn(arfidata2.*W,K,'same')./convn(W.*~isnan(arfidata),K,'same');
                    sweidata2 = sweidata;
                    sweidata2(isnan(sweidata2)) = 0;
                    W = exp(-100*(1-abs(cc)));
                    sweidata = convn(sweidata2.*W,K,'same')./convn(W.*~isnan(sweidata),K,'same');
                     
                    end
                    Cidx = 1:par.nBeams;
                    Cidx1 = Cidx;
                case 'arfiswei'
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
            
            

            if flags.calcSWS && ~strcmpi(par.mode,'arfi');
                toc1 = toc;
                fprintf('Calculating Shear Wave Velocities...')
                [TTP DX vTest] = TTPS(theta,dtheta,axial,apex,t,sweidata,sws);
                if ~iscell(sws.method);
                    swsmethod = {sws.method};
                else
                    swsmethod = sws.method;
                end
                for i = 1:length(swsmethod)
                    fprintf('%s...',sws.method{i});
                    switch upper(sws.method{i})
                        case 'TTPS'
                            CT_TTPS = sws_TTPS(TTP,DX);
                        case 'DTTPS'
                            CT_DTTPS = sws_DTTPS(TTP,DX,1);
                        case 'SMURF'
                            CT_SMURF = sws_SMURF(TTP,DX,1,par);
                    end
                end
                eval(sprintf('CT = CT_%s;',upper(swsmethod{1})));
                fprintf('done (%0.2fs)\n',toc-toc1); 
                SWS(:,:,frmii) = nanmedian(CT,3);
            end
                        
            if flags.displayflag
                toc1 = toc;
                fprintf('Generating Display...')
                realtime_display_v4
                fprintf('done (%0.2fs)\n',toc-toc1);
            end
            

            
            
            if flags.saveRes
                %outputdatadir = fullfile('C:\pjh7\SC2000\RES\',datestr(now,'yyyymmdd'));
                [pth] = fileparts(fname);
                outputdatadir = fullfile(pth,'res');
                
                % Save time stamped results file
                if par.numFrames>1
                    outputdatadir = fullfile(pth,'res',dispEst.method,sprintf('res_%s',timeStamp));
                    resfile(frmidx,:) = fullfile(outputdatadir,sprintf('res_%s_%03.0f.mat',timeStamp,frmidx));
                else
                    outputdatadir = fullfile(pth,'res',dispEst.method);
                    resfile = fullfile(outputdatadir,sprintf('res_%s.mat',timeStamp));
                end
                
                
                if ~exist(outputdatadir,'dir');mkdir(outputdatadir);end
                fprintf('Saving %s...\n',resfile(frmidx,:))
                cc = uint8(min(255,512*max(0,cc-0.5)));
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
                    'fsiq',par.fs*dispEst.interpFactor*1e6...
                    );
                btheta = Th1;
                if exist('sweidata','var') && flags.calcSWS
                              %if flags.savesws
                    save(resfile(frmidx,:),'arfidata','sweidata','cc','axial','theta','apex','t','bmodedata0','btheta','sc','scB','TTP','CT_*')
                              %else
                    %save(resfile(frmidx,:),'arfidata','sweidata','cc','axial','theta','apex','t','bmodedata0','btheta','sc','scB');    

                            %  end
                else
                save(resfile(frmidx,:), 'arfidata','cc','axial','theta','apex','t','bmodedata0','btheta','sc','scB');
                end
                toc
            end
            ARFIDATA(:,:,:,frmii) = squeeze(arfidata);
                if flags.keyboard
            keyboard
                end
        end
        ii = reshape(reshape(1:32,[],2)',[],1);
        %for i = 1:32;arfidata1(:,ii(i),:) = ARFIDATA(:,ii(i),:,i);end


%         if flags.savesws
%             switch lower(swifParams.imagingModeShort)
%                 case 'swei'
%                     swsfile = fullfile(outputdatadir,sprintf('sws_%s.mat',timeStamp));
%                     fprintf('Saving %s...\n',swsfile)
%                     dxdt = sws.dxdt;
%                     resid = sws.resid;
%                     save(swsfile,'dxdt','resid','axial','slice');
%                     toc
%             end
%         end
    end
    
    if isempty(Resfiles)
        Resfiles = resfile;
    else
        Resfiles = char(Resfiles,resfile);
    end
end
end


function [arfidata sweidata theta dtheta] = splitARFISWEI(arfidata0,par)
            arfidata = single(zeros(size(arfidata0,1),par.numBeamGroups,par.ensemble));
            sweidata = permute(reshape(arfidata0,size(arfidata0,1),par.numBeamGroups,par.nBeams,par.ensemble),[1 2 4 3]);
            [lat dlat] = genLatMatrix(par);
            dii = [-1 0 1];
            w = [0.25 0.5 0.25];
            W = repmat(w,[size(sweidata,1),1,size(sweidata,3)]);
            for i = 1:par.numBeamGroups;
                idx = max(1,min(size(sweidata,2),par.pushbeamNum(i)+dii));
                arfidata(:,par.pushbeamNum(i),:) = sum(sweidata(:,idx,:,i).*W,2)./sum(W,2);
            end
            theta = lat(1,:);
            dtheta = dlat;
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
        end
        case 'phased'
        txSpacing = (par.trackParams.paramList.txXyzGridParams.thetaX_DEG./par.trackParams.paramList.txXyzGridParams.P);
        lat = (lat.*repmat(txSpacing',[size(lat,1) 1]))';
        if isfield(par,'pushbeamNum')
            txPush = interp1(1:par.numBeamGroups,par.trackParams.paramList.txXyzGridParams.thetaX_DEG,par.pushbeamNum);
            dlat = lat - repmat(txPush(:),[1 size(lat,2)]);
        end
    end
end
