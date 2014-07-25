function prep14L5DataSWEI(resFile, parFile, leftBeams, rightBeams)
% This is a function that splits the 14L5 arfi files into L and right waves, filters and ensures proper formatting of the arfi data such that we have (axialXAlinesXtime format)
%resFile is the original resFile
%parFile is the original parFile
%leftBeams is a vector of the parts of the arfidata that constitue left of
%push track beams such as [1:7] in the current case
%rightBeams is a vector of the parts of the arfidata that constitue left of
%push track beams such as [9:15] in the current case
addpath /luscinia/ProstateStudy/abp19/14L5_proc/arfiProcCode/

res=load(resFile);
[ax Aline lateral time]= size(res.arfidata);
res.arfidata=res.arfidata(:,:,:,1:length(res.t));
% res.arfidata=res.arfidata(:,:,:,[1:3 9:length(res.t)]);  %%manual ts
% removal
% res.t=res.t([1:3 9:length(res.t)]); %%manual ts
% removal
[t, arfidata2] = filtArfiData(res.axial, res.t, res.arfidata, [50 1.8e3]);
arfidata2=reshape(arfidata2,size(arfidata2,1),rightBeams(end), [] , length(t));
axial=res.axial;
lat2=res.lat;
if(leftBeams~=0)
numBeams=length(leftBeams);
%left arfidata
res.arfidata=arfidata2(:,leftBeams,:,:);
arfidata=reshape(res.arfidata,size(res.arfidata,1),[],size(res.arfidata,4));
lat=lat2(leftBeams, :);

save('resS_22222222222222.mat', 'arfidata', 'lat', 't', 'axial')

clear arfidata
end
%right arfidata
arfidata=arfidata2(:,rightBeams,:,:);
arfidata=reshape(arfidata,size(arfidata,1),[],size(arfidata,4));
lat=lat2(rightBeams, :);

save('resS_11111111111111.mat', 'arfidata', 'lat', 't', 'axial')

clear res arfidata arfidata2 lat2 B A lat t axial
if(leftBeams~=0)
load(parFile)
nBeams=length(leftBeams);
trackParams.rxMultibeamParams.beamPatternP= trackParams.rxMultibeamParams.beamPatternP(leftBeams);

save('parS_22222222222222.mat')
end

load(parFile)
nBeams=length(rightBeams);
trackParams.rxMultibeamParams.beamPatternP= trackParams.rxMultibeamParams.beamPatternP(rightBeams);

save('parS_11111111111111.mat')


end


