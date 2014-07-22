function playARFI(res,clim,fps)
if ~exist('fps','var')
    fps = 10;
end
im = imagesc(res.lat,res.axial,res.arfidata(:,:,1),clim);
axis image;
th = title(sprintf('t = 0'));
xlabel('x (mm)');
ylabel('z (mm)');
for i = 1:length(res.t);
    tic;
    set(im,'cdata',res.arfidata(:,:,i));
    set(th,'string',sprintf('t = %0.2f',1e3*res.t(i)));
    drawnow
    while toc<(1/fps);end
end
