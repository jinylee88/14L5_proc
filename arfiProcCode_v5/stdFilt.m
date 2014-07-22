function Im = stdFilt(Im0,winsize,threshold);
Im = Im0;
for k = 1:size(Im0,3);
    for l = 1:size(Im0,4)
        im0 = Im0(:,:,k,l);
winsize = winsize - mod(winsize,2) + 1;
if length(winsize) == 1;
    winsize = winsize([1 1]);
end
impad = padarray(im0,(winsize-1)/2,nan,'both');
N = winsize(1)*winsize(2)-1;
imrep = nan([size(impad) N]);
idx = 1;
for i = -(winsize(1)-1)/2:(winsize(1)-1)/2
    for j = -(winsize(2)-1)/2:(winsize(2)-1)/2
        if i || j
           imrep(:,:,idx) = circshift(impad,[i j]);
           idx = idx+1;
        end
    end
end
sigma = nanstd(imrep((winsize(1)-1)/2+1:end-(winsize(1)-1)/2,(winsize(2)-1)/2+1:end-(winsize(2)-1)/2,:),0,3);
mu = nanmean(imrep((winsize(1)-1)/2+1:end-(winsize(1)-1)/2,(winsize(2)-1)/2+1:end-(winsize(2)-1)/2,:),3);
msk = (abs(im0-mu)./sigma) > threshold | isnan(im0);
immed = mediannan(im0,winsize);
im = im0;
im(msk) = immed(msk);
Im(:,:,k,l) = im;
    end
end