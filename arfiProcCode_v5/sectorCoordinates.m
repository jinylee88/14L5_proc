
function [X Z DX] = sectorCoordinates(theta,dtheta,r,apex)
for i = 1:size(dtheta,1);
[DTH(:,:,i) DR(:,:,i)] = meshgrid(dtheta(i,:),r);
end
DX = DR.*sind(DTH) - repmat(apex.*cosd(theta),[size(DR,1),1,size(DTH,3)]).*sind(DTH);
[TH R] =  meshgrid(theta,r);
X = R.*sind(TH) - repmat(apex.*cosd(theta),size(R,1),1).*sind(TH);
Z = R.*cosd(TH);
end
