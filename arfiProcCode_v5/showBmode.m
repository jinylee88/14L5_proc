apex = -3.54;
theta = asind(linspace(sind(-45),sind(45),size(I,2)));
axial = linspace(0,30,size(I,1));
bmodedata0 = abs(complex(I,Q));
bmodedata0 = max(0,min(1,(db(bmodedata0) - options.bmode.dynamicRange + options.bmode.offset)/options.bmode.dynamicRange));
if options.display.enable
figure(options.display.figureHandles.BMode)
clf
imsurf(theta,axial,apex,bmodedata0(:,:,[1 1 1]));
axis image
axis(options.display.axisLimits);    
end
resfile = '';
