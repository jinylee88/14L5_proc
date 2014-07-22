function process_latest(i)
set(0,'defaultfigurewindowstyle','docked')
addpath C:\pjh7\SC2000\arfiProcCode_v4\
cd V:\Program' Files'\Siemens\syngo\Bedrock\Startup\
fNames = dir('SWIF_AData*.bin')
if ~exist('i','var')
    i = 0;
end
if i == 0
    procArfi_v4(fNames(end).name);
elseif i>0
    procArfi_v4(fNames(i).name);
else
    procArfi_v4(fNames(end+i).name);
end