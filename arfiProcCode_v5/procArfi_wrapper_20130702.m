function procArfi_wrapper_20130702(fileidx)
    filepaths = {'/getlab/pjh7/SC2000/data/20130625_ExVivoAblation/pre',...
    '/getlab/pjh7/SC2000/data/20130625_ExVivoAblation/post',...
    '/getlab/pjh7/SC2000/data/20130626_ExVivoAblation/pre',...
    '/getlab/pjh7/SC2000/data/20130626_ExVivoAblation/post'};
    method = 'PesaventoFlux';
    D = {};
    for filepathidx = 1:4
    filepath = filepaths{filepathidx};
    swifFiles = dir(fullfile(filepath,'SWIF_AData*.bin'));
    for i = 1:length(swifFiles);
    D{end+1} = fullfile(filepath,swifFiles(i).name);
    end
    end
    swifFile = D{fileidx};
    procArfi_v5(swifFile);