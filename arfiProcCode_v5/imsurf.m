function imh = imsurf(th,r,apex,im,cax)
th = th(:);
r = r(:);
r = [r;r(end)+diff(r(end-1:end))];
th = [th;th(end)+diff(th(end-1:end))];
r = r-mean(diff(r))/2;
th = th-mean(diff(th))/2;
[TH R] = meshgrid(th,r);
x = R.*sind(TH) - apex*tand(TH);;
y = R.*cosd(TH);

IM = padarray(im,[1 1 0],'post');
IM(isnan(IM)) = -inf;
imh = surf(x,y,0*y,double(IM));
set(imh,'edgecolor','none');
axis ij;
axis image;
view([0 90]);
if exist('cax','var')
    caxis(cax);
end
