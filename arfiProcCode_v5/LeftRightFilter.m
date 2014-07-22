
function [vL vR] = LeftRightFilter(vel)
sz0 = size(vel);
sz = 2.^(max(0,nextpow2(size(vel))));
F = fft(fft(vel,sz(2),2),sz(3),3);
kth = linspace(-1,1,sz(2));
omega = linspace(-1,1,sz(3));
[W KTH] = meshgrid(omega,kth);
mskR = (KTH.*W<=0);
mskL = (KTH.*W>=0);
MSKR = repmat(permute(mskR,[3 1 2]),[size(vel,1),1,1,size(vel,4)]);
MSKL = repmat(permute(mskL,[3 1 2]),[size(vel,1),1,1,size(vel,4)]);
FL = 2*F.*MSKL;
FR = 2*F.*MSKR;
vR = real(ifft(ifft(FR,[],3),[],2));
vR = vR(:,1:size(vel,2),1:size(vel,3),:);
vL = real(ifft(ifft(FL,[],3),[],2));
vL = vL(:,1:size(vel,2),1:size(vel,3),:);
end


