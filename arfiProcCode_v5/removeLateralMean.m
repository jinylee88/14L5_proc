function [vel] = removeLateralMean(vel0)
vel = vel0 - repmat(mean(vel0,2),[1 size(vel0,2) 1 1]);
end