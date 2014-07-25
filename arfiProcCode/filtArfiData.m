function [t, arfidata] = filtArfiData(axial, t, arfidata, cutoffFreq, axFiltKer, stdCutoff);
% function [t arfidata] = filtArfiData(axial, t, arfidata, cutoffFreq, axFiltKer, stdCutoff);
%
% this function is designed to work with a single acquisition (numAcq=1)
% if multiple acquisitions were performed, feed them in one at a time
%
% if nargin<4,cutoffFreq = [20 1000];end % filter cutoff frequencies (Hz)
% if nargin<5,axFiltKer = 0.5;end % axial filter length (mm)
% if nargin<6,stdCutoff = 10;end % cutoff in microns for the minimum standard deviation for push/reverb time steps

if nargin<4,cutoffFreq = [20 1000];end % filter cutoff frequencies (Hz)
if nargin<5,axFiltKer = 0.5;end % axial filter length (mm)
if nargin<6,stdCutoff = 30;end % cutoff in microns for the minimum standard deviation for push/reverb time steps


% determine push and reverb time steps, then remove them
if ndims(arfidata)==3
    arfidata=arfidata(:,:,1:length(t));
    out = squeeze(max(std(arfidata,0,2),[],1));
    ts = out<stdCutoff; % valid time steps
    if sum(ts)<10 || sum(ts)>(length(out)-4)
        out2=sort(-out);
        cut=-out2(4);
        ts=out<cut;
        ts(1)=1;
        outEnd=length(out);
        ts(outEnd)=1;
    end
elseif ndims(arfidata)==4
    arfidata=arfidata(:,:,:,1:length(t));
    out = squeeze(max(max(std(arfidata,0,3),[],1),[],2));
    ts = out<stdCutoff; % valid time steps
    if sum(ts)<10 || sum(ts)>(length(out)-4)
        out2=sort(-out);
        cut=-out2(4);
        ts=out<cut;
        ts(1)=1;
         outEnd=length(out);
        ts(outEnd)=1;
    end
end

% Set up filtering parameters
dt = min(diff(t));
if (t(end)-t(end-1))>10*dt,ts(end)=0;end
fprintf(1, 'Removing time step:\t%d\n', find(ts==0));
t = t(ts);t = round(t*1e4)/1e4; % numerical tolerance issues
tn = t(1):dt:t(end);tn = round(tn*1e4)/1e4; % numerical tolerance issues
fs = 1./dt*1e3;
[B A] = butter(2, cutoffFreq./(fs/2)); % filter at 20Hz-1kHz
n = max(1,round(axFiltKer./mean(diff(axial)))); % axial filter (minimum 1 sample)

% interpolate and filter the data
arfidata = reshape(arfidata, size(arfidata,1), [], size(arfidata,ndims(arfidata)));
arfidata = arfidata(:,:,ts);

tstart = tic;
arfidata = medfilt1(double(arfidata), double(n), [], 1); % axial filter
fprintf(1, 'Axial median filter complete in %0.2f seconds\n', toc(tstart));

tstart = tic;
arfidata = temporalFilter(arfidata, t, tn, B, A); % interpolate and filter in time
fprintf(1, 'Temporal filter complete in %0.2f seconds\n', toc(tstart));

t = tn;

end

function data = temporalFilter(arfidata, t, tn, B, A)
% interpolate and filter
D = size(arfidata);
arfidata = reshape(arfidata, [], D(end));
data = interp1(t(:), arfidata', tn(:), 'spline')';
D(end) = length(tn);
data = single(filtfilt(B,A,data'))';
data = reshape(data,D);
end
