function cmap = Hollender(N)
if ~exist('N','var')
    N = 64;
end

up = linspace(0,1,floor(N/4)+1);
up = up(2:end)';
dn = linspace(1,0,floor(N/4)+1);
dn = dn(1:end-1)';
z = zeros(floor(N/4),1);
o = ones(floor(N/4),1);

up1 = linspace(0,1,round(N/4)+1);
up1 = up1(2:end)';
dn1 = linspace(1,0,round(N/4)+1);
dn1 = dn1(1:end-1)';
z1 = zeros(round(N/4),1);
o1 = ones(round(N/4),1);

if mod(N,2) == 1
cmap = [o1 dn1 z1;dn z z;0 0 0;z z up; z1 up1 o1;];
else
cmap = [o1 dn1 z1;dn z z;z z up; z1 up1 o1;];
end   

cmap = flipud(cmap);
trim = size(cmap,1)-N;
cmap = cmap(floor(trim/2)+(1:N),:);