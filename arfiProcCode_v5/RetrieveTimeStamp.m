function timestamp = RetrieveTimeStamp(filename)
if isempty(filename)
    timestamp = '';
else
for i = 1:size(filename,1)
timeidx = regexp(filename(i,:),'\d\d\d\d\d\d\d\d\d\d\d\d\d');
endidx = timeidx +13;
timestamp(i,:) = char(filename(i,timeidx:endidx));
end
end