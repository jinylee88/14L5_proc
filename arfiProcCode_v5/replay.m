figure(1);
clf
im = imsurf(theta,axial,apex,sweidata(:,:,1,1),[-10 10]);
colormap(Hollender)
for j = 1:32;
    for i = 1:size(sweidata,3)-1;
        set(im,'cdata',padarray(double(sweidata(:,:,i,j)),[1 1],'post'));
        pause(0.01);
    end;
end
