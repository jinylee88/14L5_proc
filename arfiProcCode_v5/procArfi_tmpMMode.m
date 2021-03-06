function Resfiles = procArfi_v4(varargin)
tic

flags = struct(...
    'calcSWS',         1 ...
    ,'swsdisplayflag',  0 ...
    ,'saveRes',        1 ...
    ,'displayflag',     0 ...
    ,'savesws',         1 ...
    ,'cc',              1 ...
    ,'keyboard',         0 ...
    );

mf = struct(...
    'Enable',       1,...
    'Order',        2,...
    'tMotion_ms',   [-inf -0.15 2.0 inf]);

rt = struct(...
    'dispLim',     [0 20],...
    'ccRng',    [0.95 1],...
    'axisLimits',    [-7.5 7.5 0 25],...
    't_arfi_ms', 0.6, ...
    'bdB',      30,...
    'bdBoffset',0,...
    'bMask',1,...
    'FrameIndices','all',...
    'depthCorrect',0);

dispEst = struct(...
    'method','kasai',...
    'interpFactor',1,...
    'kernelLength',3,...
    'searchRegion',1,...
    'removeKickBack',0);

sws = struct(...
    'dxdtrange',    [0.5 8]...
    ,'xrange',       [1 5]...
    ,'trange',       [0 inf]...
    ,'smoothingkernel',[5 1]);
sws.method = {'DTTPS','TTPS'};

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
        toc1 = toc;
        fprintf('Setting up Paths...')
        s = path;
        if ~any(strfind(s,'C:\pjh7\SC2000\arfiProcCode_v4'));
            addpath C:\pjh7\SC2000\arfiProcCode_v4
        end
        if ~any(strfind(s,'C:\pjh7\SC2000\scanconversion'));
            addpath C:\pjh7\SC2000\scanconversion\
        end
        set(0,'defaultaxesfontweight','bold')
        [datadir] = fileparts(fname);
        if isempty(datadir)
            datadir = pwd;
        end
        tmpdir = 'C:\pjh7\SC2000\tmp\';
        fprintf(' done (%0.2fs)\n',toc-toc1);
    end
    
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
            priUsec = [max(par.priusec)*ones(1,par.nref+length(par.pushFocalDepth))];
            for i = 1:length(par.priusec);
                priUsec = [priUsec par.priusec(i)*ones(1,par.ntrack(i))];
            end
            priUsec = [-1*sum(priUsec(1:par.nref+length(par.push_focaldepth)-1)) priUsec(1:end-1)];
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
        clear Bmodedata0 Bimg ARFIDATA cdata SWS
        
        
        fs = par.fs*1e6;
        fc = par.trackParams.reqList.txPulseReqs.TxPulse.fm_mhz*1e6; % Hz
        c = par.c; % m/s
        NumIter = 3;
        
        I0 = I(:,:,1);
        Q0 = Q(:,:,1);
        D = size(I0);
        D(1) = D(1).*dispEst.interpFactor;
        if dispEst.interpFactor>1
            [Iup0,Qup0] = computeUpsampledIQdata(I0,Q0,dispEst.interpFactor);
            Iup0 = reshape(Iup0, D);
            Qup0 = reshape(Qup0, D);
        else
            Iup0 = I0;
            Qup0 = Q0;
        end
        fs = fs*dispEst.interpFactor;
        KLen = ceil(dispEst.kernelLength*fs/fc);
        if KLen < 3
            warning('Extending Kernel Length to 3 sample (%0.1f wavelengths).',KLen*fc/fs);
            KLen = 3;
        end
        SrchLen = round(dispEst.searchRegion*fs/fc);
        fdem = par.trackParams.Apl3.Mod(1).DsF.data*1e6; % frequency dataset values (MHz shift)
        frange = par.trackParams.Apl3.Mod(1).DsF.rr;  % reference ranges (mm)
        fc_vec = reshape(interp1(frange, fdem, axial0), size(Iup0,1), 1);
        fdem_vec = fc_vec./fs;
        
        IQ0 = complex(Iup0,Qup0);
        
        nFrames = size(I,3);
        u = single(zeros(size(IQ0,1),size(IQ0,2),nFrames-1));
        cc = u;
        b = u;
        uk = zeros(1,nFrames-1);
        cck = zeros(1,nFrames-1);
        
        tstart = tic;
        
        fprintf('Calculating Displacements...\n')
                
        %s = sprintf('%5.0f/%5.0f (Est. Time Remaining %3.0f:%02.0f)',0,nFrames-1,0,0);
        %fprintf('%s',s);
        interpFactor = dispEst.interpFactor;
        I00 = I(:,:,1:nFrames-1);
        Q00 = Q(:,:,1:nFrames-1);
        I01 = I(:,:,2:nFrames);
        Q01 = Q(:,:,2:nFrames);
        
        parfor refidx = 1:nFrames-1;
%             elapsedTime = toc(tstart);
%             estTimeRemaining = (elapsedTime/refidx)*(nFrames-refidx);
%             Tsec = mod(estTimeRemaining,60);
%             Tmin = (estTimeRemaining-Tsec)/60;
%             fprintf('%s',8*ones(1,length(s)));
%             s = sprintf('%5.0f/%5.0f (Est. Time Remaining %3.0f:%02.0f)',refidx,nFrames-1,Tmin,floor(Tsec));
%             fprintf('%s',s);



            %fprintf('%s%5.0f/%5.0f',8*ones(1,11),refidx,nFrames-1);
            
              %I0 = I(:,:,refidx);
              %Q0 = Q(:,:,refidx);
              I0 = I00(:,:,refidx);
              Q0 = Q00(:,:,refidx);
              
        if interpFactor>1
            [Iup0,Qup0] = computeUpsampledIQdata(I0,Q0,interpFactor);
            Iup0 = reshape(Iup0, D);
            Qup0 = reshape(Qup0, D);
        else
            Iup0 = I0;
            Qup0 = Q0;
        end
            IQ0 = complex(Iup0,Qup0);
            
            %I1 = I(:,:,refidx+1);
            %Q1 = Q(:,:,refidx+1);
            I1 = I01(:,:,refidx);
            Q1 = Q01(:,:,refidx);
            
            if interpFactor>1
                [Iup,Qup] = computeUpsampledIQdata(I1,Q1,interpFactor);
                Iup = reshape(Iup, D);
                Qup = reshape(Qup, D);
            else
                Iup = I1;
                Qup = Q1;
            end
            
            IQ = complex(Iup,Qup);
            
            [Tau cc0]= PesaventoParallel3(IQ0,IQ,fs,fc,KLen,SrchLen,NumIter);
            
            u(:,:,refidx) = Tau*c/2*1e6;
            cc(:,:,refidx) = cc0;
            b(:,:,refidx) = db(abs(IQ0));
            
            %if dispEst.removeKickBack
                [TauK ccK]= PesaventoParallelWholeLine2(IQ0,IQ,fs,fc,SrchLen,NumIter);
                uk(refidx) = TauK*c/2*1e6;
                cck(refidx) = ccK;
            %end
            
            %IQ0 = IQ;
        end
        fprintf(1, '\nDisplacement Computation Time: %0.2fs\n', toc(tstart));
        
       
        t = (0:nFrames-1)*par.priusec*1e-6;
        r = axial0;
        apex = par.trackParams.specList.scanSpecs.apexMm(3);
        theta = genLatMatrix(par);
        par.imagingdepth = 0.1*max(axial0);
        lambdamicron = ((c/2)/fc)*1e6;
        uPhaseUnwrap = u;
        uPhaseUnwrap(u>lambdamicron/2) = uPhaseUnwrap(u>lambdamicron/2)-lambdamicron;
        uPhaseUnwrap(u<-lambdamicron/2) = uPhaseUnwrap(u<-lambdamicron/2)+lambdamicron;
        
        if flags.saveRes
            [pth] = fileparts(fname);
            outputdatadir = fullfile(pth,'res');
            resfile = fullfile(outputdatadir,sprintf('u_%s.mat',timeStamp));
            
            if ~exist(outputdatadir,'dir');mkdir(outputdatadir);end
            fprintf('Saving %s...\n',resfile)
            
            sc = struct(...
                'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
                'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
                'latinc',.1e-3,...
                'axialmin',-1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
                'axialmax',1e-2*(par.imagingdepth)- 1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
                'axialinc',.1e-3,...
                'min_phi',min(theta),...
                'span_phi',max(theta)-min(theta),...
                'apex',0.1*par.trackParams.paramList.txGridParams.apexMm(3),...
                'fsiq',par.fs*dispEst.interpFactor*1e6...
                );
             cc = uint8(-100*db(abs(cc)));  
             b = uint8(255*min(1,max(0,(b+rt.bdBoffset-rt.bdB)/rt.bdB)));
             u = int16(2^16*(uPhaseUnwrap/lambdamicron));
             %if dispEst.removeKickBack
                 cck = uint8(-100*db(abs(cck)));
                  ukPhaseUnwrap = uk;
                  ukPhaseUnwrap(uk>lambdamicron/2) = ukPhaseUnwrap(uk>lambdamicron/2)-lambdamicron;
                  ukPhaseUnwrap(uk<-lambdamicron/2) = ukPhaseUnwrap(uk<-lambdamicron/2)+lambdamicron;
                  uk = int16(2^16*(ukPhaseUnwrap/lambdamicron));
             %end

                 
            save(resfile,'u','b','cc','t','r','theta','uk','cck','sc','lambdamicron','apex');
            
            toc
        end
        if flags.keyboard
            keyboard
        end
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
