function phantomPar = readCompressARFIDyn(filename);
fid = fopen(filename,'r');
s = fgetl(fid);
s = fgetl(fid);
C = textscan(fid,'%s%s%s',1);
phantomPar.Geometry = C{3}{1};
while ~strcmpi(s,'$$$$ BACKGROUND MATERIAL $$$$')
    s = fgetl(fid);
end
s = fgetl(fid);
C = textscan(fid,'%f%f%f%f%f%f%f',1,'delimiter',',');
phantomPar.BGkPa = C{3}/1e4;
phantomPar.BGpoissons = C{4};
while ~strcmpi(s,'$$$$ SPHERE MATERIAL $$$$')
    s = fgetl(fid);
end
s = fgetl(fid);
C = textscan(fid,'%f%f%f%f%f%f%f',1,'delimiter',',');
phantomPar.FGkPa = C{3}/1e4;
phantomPar.FGpoissons = C{4};
fclose(fid);