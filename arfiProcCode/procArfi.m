function procArfi(fname, parFile, interpFactor, kernelLength, ccmode, ref_idx)
% function procArfi(fname, parFile, interpFactor, kernelLength, ccmode)
%
% Inputs: fname - file name of binary file saved using analytic pipeline (SWIF_AData*.bin)
%         parFile - parameters file saved when arfi_image.m was run (parameters.mat)
%         interpFactor - upsampling factor to use (5)
%         kernelLength - length of kernel to use in loupas in wavelengths (4)
%         ccmode - whether or not to compute correlation coeffients (0)
%         ref_idx - index of time step to use as the reference for displacement estimation (par.nref)
% All inputs are optional and will default to the values in parentheses if ommitted

addpath(fileparts(which(mfilename)));

% Check inputs and set default parameters
if nargin<1
    fname = dir('SWIF_AData*.bin');fname = fname(end).name;
end
% Pull out time stamp from filename
[basePath,fname] = fileparts(fname);
timeStamp = fname(12:end);
fname = fullfile(basePath, strcat(fname, '.bin'));
if nargin<2 || isempty(parFile)
    parFile = fullfile(pwd, sprintf('par_%s.mat',timeStamp));
    % Check to see if there is a time stamped parameters file, if it is in current directory, or one directory level higher
    if ~exist(parFile, 'file'), parFile = fullfile(pwd, 'parameters.mat');end
    if ~exist(parFile, 'file'), parFile = fullfile(pwd, '..', 'parameters.mat');end 
end
if nargin<3
    interpFactor = 5;
end
if nargin<4
    kernelLength = 4;
end
if nargin<5
    ccmode = 0;
end
if nargin<6
    ref_idx = []; % unless explicitly supplied, value of ref_idx is set after loading in par file and detecting location of push frame
end

% Check to make sure data, dimensions, and parameters files all exist
dimsname = ['SWIF_ADataDims_' timeStamp '.txt'];
if ~exist(fname, 'file'),error('No data file found');end
if ~exist(dimsname, 'file'),error('No dimensions file found');end
if ~exist(parFile, 'file'), error('No parameters file found');end

% Pull out IQ data
data = readSwif(fname, dimsname);
I = single(data.I);
Q = single(data.Q);
clear data

% Load parameters file
par = load(parFile);
par.kernelLength = kernelLength;
par.interpFactor = interpFactor;

% Confirm location of push frame in case of DMA-SWIF Buffer Event
temp_cc = computeCC(I(:,round(size(I,2)/2),1:par.nref+par.npush),round(size(I,1)/2));
temp_cc = squeeze(mean(temp_cc,1));
temp_cc(find(isnan(temp_cc))) = 0;
acq_nref = find(abs(temp_cc)<0.5,1)-1;
if acq_nref == par.nref
    fprintf(1,'Push frame at expected location\n');
elseif acq_nref == par.nref-1
    fprintf(1,'Push frame offset due to DMA-SWIF Buffer Event \npar.nref, par.ensemble will be reduced by 1\n');
    par.nref = par.nref-1; par.ensemble = par.ensemble-1;
    I = I(:,:,1:end-1); Q = Q(:,:,1:end-1);
else
    warning('Unable to detect location of push frame using correlation coefficients. Check IQ Data')
end

if isempty(ref_idx)
    ref_idx = par.nref;
end
par.ref_idx = ref_idx;

% add number of acquisitions and harmonic flag to parameters structure if it isn't there
if ~isfield(par, 'numAcq'),par.numAcq = 1;end
if ~isfield(par, 'isHarmonic'), par.isHarmonic = 0;end

% % Reshape and reorder the data if we are using multiple acquisitions
if par.numAcq>1
    I=reshape(I,size(I,1),size(I,2),size(I,3)/par.numAcq,par.numAcq);
    I=reshape(permute(I,[1,2,4,3]),size(I,1),size(I,2)*par.numAcq,size(I,3));
    Q=reshape(Q,size(Q,1),size(Q,2),size(Q,3)/par.numAcq,par.numAcq);
    Q=reshape(permute(Q,[1,2,4,3]),size(Q,1),size(Q,2)*par.numAcq,size(Q,3));
end

% sum data as needed if we are doing harmonic acquisitions
if par.isHarmonic
    [I,Q] = genHarmonicSummedData(I,Q,par);
end

% Compute axial vector
N = size(I,1)*interpFactor;
axial = (0:N-1)*(par.c/1e3)/(2*par.fs*interpFactor);

% find center frequency (double frequency for harmonic data)
if par.isHarmonic
    par.fc = par.trackParams.fc*1e6*2; % Hz
else
    par.fc = par.trackParams.fc*1e6; % Hz
end
par.lambda = par.c / par.fc * 1e3; % mm

% Compute displacements using the last reference and then reorder the data
if par.ref_idx == -1
    fprintf(1,'Computing displacements: Progressive\n')
else
    fprintf(1,'Computing displacements: Anchored at Frame %d (par.nref = %d)\n',par.ref_idx,par.nref)
end

[arfidata, I, Q] = runLoupas(I, Q, interpFactor, kernelLength, axial, par);
arfidata = single(arfidata);

if ccmode
    fprintf(1, 'Computing complex correlation coefficients\n');
    fs = par.fs*1e6*interpFactor;
    fc = par.fc; % Hz
    cc_coef = single(abs(computeCC(complex(I,Q),round(kernelLength*fs/fc),1)));
end

% Remove extra axial samples
axial = axial(1:size(arfidata,1));
axial = axial + par.kernelLength * par.lambda / 2; % shift axial vector based on 50% of tracking kernel

% generate time vector
% time = 0 is the first tracking vector
[t, txTypeIndex, pushPRF] = genTimeVector(par);
par.txTypeIndex = txTypeIndex;
par.pushPRF = pushPRF;

% generate lateral position vector
lat = genLatMatrix(par);
if par.separateFocalZoneAcqs
    arfidata = reshape(arfidata, [size(arfidata,1), par.nBeams, par.numBeamGroups*length(par.pushFocalDepth), par.numAcq, size(arfidata,3)]); % reshape into 4-D matrix
    if ccmode,cc_coef = reshape(cc_coef, [size(cc_coef,1), par.nBeams, par.numBeamGroups*length(par.pushFocalDepth), par.numAcq, size(cc_coef,3)]);end
else
    arfidata = reshape(arfidata, [size(arfidata,1), par.nBeams, par.numBeamGroups, par.numAcq, size(arfidata,3)]); % reshape into 4-D matrix
    if ccmode,cc_coef = reshape(cc_coef, [size(cc_coef,1), par.nBeams, par.numBeamGroups, par.numAcq, size(cc_coef,3)]);end
end

switch lower(par.mode)
    case {'mmode'}
        arfidata(:,:,par.trackBeamNum,:,:) = arfidata;
        if ccmode,cc_coef(:,:,par.trackBeamNum,:,:) = cc_coef;end
end

% collapse to 3D matrix if lat variable is a vector
if sum(size(lat)~=1)==1
    arfidata = reshape(arfidata, size(arfidata,1), [], size(arfidata,4), size(arfidata,5));
    if ccmode,cc_coef = reshape(cc_coef, size(cc_coef,1), [], size(cc_coef,4), size(cc_coef,5));end
    if par.numAcq>1
        arfidata = reshape(permute(arfidata, [1 2 4 3]), size(arfidata,1), size(arfidata,2), size(arfidata,4), 1, size(arfidata,3));
        if ccmode,cc_coef = reshape(permute(cc_coef, [1 2 4 3]), size(cc_coef,1), size(cc_coef,2), size(cc_coef,4), 1, size(cc_coef,3));end
    end
elseif par.numAcq>1
    arfidata = permute(arfidata, [1 2 3 5 4]);
    if ccmode,cc_coef = permute(cc_coef, [1 2 3 5 4]);end
end
arfidata = squeeze(arfidata);
if ccmode,cc_coef = squeeze(cc_coef);end

% read in b-mode data if it was acquired
if par.BmodeAcq
    bmodeParFname = sprintf('SWIF_BModeOutputImageDims0_%s.txt', timeStamp);
    bmodeFname = sprintf('SWIF_BModeOutputImage0_%s.img', timeStamp);
    if ~exist(bmodeParFname, 'file')
        warning('B-mode parameters file not detected, B-mode images may not have been saved');
        bmodeSave = 0;
    elseif ~exist(bmodeFname, 'file')
        warning('B-mode data file not detected, B-mode images may not have been saved');
        bmodeSave = 0;
    else
        bmodePar = struct;
        bmodeParFid = fopen(bmodeParFname, 'r');
        while 1
            tline = fgetl(bmodeParFid);
            if ~ischar(tline)||isempty(tline), break, end
            k = strfind(tline, ':');
            if ~isempty(k)
                eval(sprintf('bmodePar.%s = %s;', tline(1:k-1), tline(k+1:end)));
            else
                if isempty(which('strsplit'))
                    tmp = textscan(tline,'%s');
                    tmp = tmp{1};
                    try
                        eval(sprintf('bmodePar.%s = %s;', tmp{1}, tmp{2}));
                    catch err
                        try
                            eval(sprintf('bmodePar.%s = ''%s'';', tmp{1}, tmp{2}));
                        catch
                            if length(tmp)==2
                                rethrow(err)
                            end
                        end
                    end
                else
                    tmp = strsplit(tline);
                    try
                        eval(sprintf('bmodePar.%s = %s;', tmp{2}, tmp{3}));
                    catch err
                        try
                            eval(sprintf('bmodePar.%s = ''%s'';', tmp{2}, tmp{3}));
                        catch
                            rethrow(err)
                        end
                    end
                end
            end
        end
        fclose(bmodeParFid);
        bmodeFid = fopen(bmodeFname, 'rb');
        bimg = fread(bmodeFid, inf, '*uint8');
        fclose(bmodeFid);
        bimg = reshape(bimg, bmodePar.SamplesPerLine, bmodePar.LinesPerSlice, []);
        bax = (0:bmodePar.SamplesPerLine-1)./bmodePar.NumSamplesPerMm;
        blat = linspace(bmodePar.FirstLinePos, bmodePar.LastLinePos, bmodePar.LinesPerSlice);
        bmodeSave = 1;
    end
else
    bmodeSave = 0;
end

% change axial and lat to radial and angular for phased and curvilinear
% array
if max(strcmp(par.probeType,{'curvilinear','phased'}))
    radial = axial;
    angular = lat;
    apex = par.trackParams.txXyzGridParams.apexMm(3);
    clear axial lat 
end

% Check existing parameters, add proc parameters to parFile
par = checkParams(par);
save(parFile,'-struct','par');

% Save time stamped results file
resfile = ['res_' timeStamp '.mat'];
if strcmp(par.probeType,'linear')
    if bmodeSave
        if exist('cc_coef', 'var')
            save(resfile, 'arfidata', 'axial', 'lat', 't', 'cc_coef', 'bimg', 'blat', 'bax', '-v7.3');
        else
            save(resfile, 'arfidata', 'axial', 'lat', 't', 'bimg', 'blat', 'bax', '-v7.3');
        end
    else
        if exist('cc_coef', 'var')
            save(resfile, 'arfidata', 'axial', 'lat', 't', 'cc_coef', '-v7.3');
        else
            save(resfile, 'arfidata', 'axial', 'lat', 't', '-v7.3');
        end
    end
elseif max(strcmp(par.probeType,{'curvilinear','phased'}))
    if bmodeSave
        if exist('cc_coef', 'var')
            save(resfile, 'arfidata', 'radial', 'angular', 'apex', 't', 'cc_coef', 'bimg', 'blat', 'bax', '-v7.3');
        else
            save(resfile, 'arfidata', 'radial', 'angular', 'apex', 't', 'bimg', 'blat', 'bax', '-v7.3');
        end
    else
        if exist('cc_coef', 'var')
            save(resfile, 'arfidata', 'radial', 'angular', 'apex', 't', 'cc_coef', '-v7.3');
        else
            save(resfile, 'arfidata', 'radial', 'angular', 'apex', 't', '-v7.3');
        end
    end
else
    error('Unknown probe type')
end
    
    
    
    
    
