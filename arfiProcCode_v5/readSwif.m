function [swif swifParams] = readSwif(fname, dimsname);
% function [swif swifParams] = readSwif(fname, dimsname);
%
% Function to read in binary data saved using the analytic pipeline on the SC2000
% 
% Inputs: fname - file name of data (optional)
%         dimsname - file name of dimensions file (optional)
% 
% Will automatically use the last data and dimensions files based on name if no inputs are given
%
% Outputs: swif - structure that contains 2 fields: I and Q
%          swifParams - structure that contains all of the parameters read in from the dimensions file
%
% sjr6 3/12/12
% original version: sjr6 while interning at Siemens (7/1/11)

if nargin<1
    fname = dir('SWIF_AData*.bin');fname = fname(end).name;
    if ~exist(fname, 'file'),error('No data file found');end
end
if nargin<2
    dimsname = dir('SWIF_ADataDims_*.txt');dimsname = dimsname(end).name;
    if ~exist(dimsname, 'file'),error('No dimensions file found');end
end

fprintf(1, 'filename: %s\n', fname);
fid = fopen(dimsname, 'r');
swifParams = struct;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        param = textscan(tline, '%s'); param = param{1};

        isHex = false;
        for i = 1:length(param)
            if (strcmpi(param(i), 'pointer'))
                isHex = true;
                break
            end
        end

        if (~isHex)
            val = str2mat(param{end});
            if strcmp(param{1},'focus')
                val = ['[' strrep(strrep(param{end}(2:end-2),'p','.'),'_',' ') ']'];
            end
        else
            val = mat2str(hex2dec(param{end}));
        end

        param = [param{1:end-2}];
        if (~isfield(swifParams, param))
            eval(sprintf('swifParams.%s = %s;', param, val));
        elseif (isfield(swifParams, param) && isempty(eval(sprintf('swifParams(end).%s', param))))
            eval(sprintf('swifParams(end).%s = %s;', param, val))
        else
            eval(sprintf('swifParams(end+1).%s = %s;', param, val));
        end
    end
end
fclose(fid);

fid = fopen(fname, 'rb');
data = fread(fid, inf, 'int16');
fclose(fid);

data = reshape(data, swifParams(1).SamplesPerLine*2, [], swifParams.FrameCount);
swif.I = data(1:2:end, :, :);
swif.Q = data(2:2:end, :, :);