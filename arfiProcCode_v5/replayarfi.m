figure(1);
clf
im = imsurf(theta,axial,apex,arfidata(:,:,1),[0 10]);
colormap(jet)
    for i = 1:size(arfidata,3);
        set(im,'cdata',padarray(double(arfidata(:,:,i)),[1 1],'post'));
        pause(0.1);
    end;
