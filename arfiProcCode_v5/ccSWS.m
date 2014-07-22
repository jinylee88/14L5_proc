function sws = ccSWS(t,x,slices)
slice = diff(slices(:,:,1),1,2)/mean(diff(t));
ctest = 0.5:0.025:15;
v = mean(slice(floor((size(slice,1)+1)/2):ceil((size(slice,1)+1)/2),:),1);
t0 = t(1:end-1);
for i = 1:length(ctest)
    dt = (abs(x)/ctest(i));
    for j = 1:length(x);
        tst(j,:,i) = interp1(t0,v,t0-dt(j));
    end
end
for j = 1:size(slices,3);
    slice = diff(slices(:,:,j),1,2)/mean(diff(t));
    for i = 1:length(ctest)
        ncc(j,i) = nanmean(nanmean(tst(:,:,i).*slice));
    end
end
%[pk idx] = max(ncc,[],2);
%sws = ctest(idx);
[pk sws] = subsamplepeak(ctest,ncc,2);
sws = sws';