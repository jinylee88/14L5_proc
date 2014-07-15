function svm=mom1v(vec,length,window,svm)

%
% routine to calculate the mean within a sliding window
%

% This function was copied from:
% /nefs/ncr2/ProcessArfiData/v3.0/routines/mom1v.m
%
% sjr6 3/12/12

svm(1)=0;           % get the first sum over window
for i=1:window
    svm(1) = svm(1)+vec(i);
end
                        % for each subsequent window,
                        %   subtract first value and add next value
for i=window+1:length
    svm(i-window+1) = svm(i-window)-vec(i-window)+vec(i);
end
svm=svm/window;     % get mean over window size
