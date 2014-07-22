function paths = getSWIFdir(rootpath);
swifFiles = dir(fullfile(rootpath,'SWIF_AData*.bin'));
if ~isempty(swifFiles);
    paths = rootpath;
else
    d = dir(rootpath);
    if length(d)<3
        paths = '';
    else
        
        for i = 3:length(d);
            pathi = getSWIFdir(fullfile(rootpath,d(i).name));
            if ~isempty(pathi)
            if exist('paths','var')
            paths = char(paths,pathi);
            else
            paths = pathi;
            end
            end
        end
        if ~exist('paths','var')
            paths = '';
        end
    end
end