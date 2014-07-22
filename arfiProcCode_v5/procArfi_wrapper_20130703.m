function procArfi_wrapper_20130703(fileidx)
    filepaths = {...
    '/getlab/pjh7/SC2000/data/20130219_RA_PostAblation/',...
    '/getlab/pjh7/SC2000/data/20130219_RA_PreAblation/'...
    };
    method = 'PesaventoFlux';
    D = {};
    for filepathidx = 1:2
    filepath = filepaths{filepathidx};
    swifFiles = dir(fullfile(filepath,'SWIF_AData*.bin'));
    for i = 1:length(swifFiles);
    D{end+1} = fullfile(filepath,swifFiles(i).name);
    end
    end
    swifFile = D{fileidx};
    procArfi_v5(swifFile);
    
    