function [CT TTP] = MModeTTPS(uData,swsOptions);
theta = uData.theta(1,:);
dtheta = theta;
axial = uData.r;
apex = uData.apex;
t1 = uData.t(1:55)*1e3;
[X Z DX] = sectorCoordinates(theta,dtheta,axial,apex);
DX3 = repmat(permute(DX,[1 2 4 3]),[1 1 length(t1) 1]);
T3 = repmat(reshape(t1(:),1,1,[]),[size(DX3,1),size(DX3,2),1,size(DX3,4)]);
DXDT = abs(DX3)./T3;
nanMatrix = ones(size(uData.u,1),size(uData.u,2),55);
nanMatrix(DXDT<swsOptions.dxdtrange(1)) = nan;
nanMatrix(DXDT>swsOptions.dxdtrange(2)) = nan;
nanMatrix(abs(DX3)<swsOptions.xrange(1)) = nan;
nanMatrix(abs(DX3)>swsOptions.xrange(2)) = nan;
nanMatrix(T3<swsOptions.trange(1)) = nan;
nanMatrix(T3>swsOptions.trange(2)) = nan;
fprintf('Frame %3.0f/192',0)
TTP = zeros(length(axial),length(theta),192);
CT =  zeros(length(axial),length(theta),192);
for frameIdx = 1:192
    fprintf('%s%3.0f/192',8*ones(1,7),frameIdx);
    vel = double(uData.u(:,:,(frameIdx-1)*56 + (1:55)))*uData.lambdamicron/(2^16);
    cc = ones(size(vel));
    cc(uData.cc(:,:,(frameIdx-1)*56 + (1:55))>1) = nan;
    [pk TTP(:,:,frameIdx)] = subsamplepeak(t1,vel.*nanMatrix.*cc,3);
    CT(:,:,frameIdx) = abs(DX)./TTP(:,:,frameIdx);
end
fprintf('...done\n')
