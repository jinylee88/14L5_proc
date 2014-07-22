    good = false(1,length(D));
for i = 1:length(D)
    swifFile = D{i};
    filepath = fileparts(swifFile);
    timestamp = RetrieveTimeStamp(swifFile);
    respath = fullfile(filepath,'res','PesaventoFlux');
    resFile = fullfile(respath,sprintf('res_%s.mat',timestamp));
    good(i) = exist(resFile,'file');
end
%%
badidx = find(~good);
for i = 1:length(badidx)
    swifFile = D{badidx(i)};
    [pth fname ext] = fileparts(swifFile);
    timestamp = RetrieveTimeStamp(fname);
    dimsfile = sprintf('SWIF_ADataDims_%s.txt',timestamp);
    d0 = dir(fullfile(pth,dimsfile));
    d = dir(fullfile(pth,'SWIF_ADataDims*.txt'));
    filesize = median([d.bytes]);
    okidx = find([d.bytes]==filesize);
    if isempty(d0) || d0.bytes<filesize;
        copyfile(fullfile(pth,d(okidx(1)).name),fullfile(pth,'tmp.txt'),'f');
        movefile(fullfile(pth,'tmp.txt'),fullfile(pth,dimsfile));
    end
    procArfi_v5(D{badidx(i)});
end
