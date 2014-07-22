function [vel] = removeLateralMedian(vel0)
vel = vel0 - repmat(median(vel0,2),[1 size(vel0,2) 1 1]);
end