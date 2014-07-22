function [pk tpk A]=subsamplepeak(t,data,dim)
t = t(:);
if ~exist('dim','var')
    dim = 1;
end
[pk0 idx0] = nanmax(data,[],dim);
t0 = reshape(t(idx0),size(idx0));
t00 = reshape(t(max(idx0-1,1)),size(idx0));
t01 = reshape(t(min(idx0+1,length(t))),size(idx0));
siz = size(data);
nd = length(siz);
for i = 1:nd
    a{i} = 1:siz(i);
end
outstr = sprintf('%s',eval(sprintf('[%s ]',sprintf('char(64+%g),'' '' ',[1:nd]))));
cmdstr = sprintf('%s',sprintf('[%s\b] = ndgrid(%s\b);',outstr,sprintf('a{%g},',[1:nd])));
bkspc = (double(cmdstr) == 8);
bkspc(find(bkspc)-1) = 1;
cmdstr = cmdstr(~bkspc);
eval(cmdstr);
pk00 = pk0;
pk01 = pk0;
AA = cell(1,nd);
for i = 1:nd;
    eval(sprintf('AA{%g} = mean(%s,dim);',i,char(64+i)));
end
AA{dim} = max(idx0-1,1);
varstr = sprintf('AA{%g}(:),',[1:nd]);
varstr = varstr(1:end-1);
cmdstr = sprintf('IND00 = sub2ind(siz,%s);',varstr);
eval(cmdstr);
AA{dim} = min(length(t),idx0+1);
cmdstr = sprintf('IND01 = sub2ind(siz,%s);',varstr);
eval(cmdstr);
pk00(:) = data(IND00);
pk01(:) = data(IND01);


  d = (t00 - t0) .* (t00 - t01) .* (t0 - t01);%denominator
  a = (t01  .* (pk0 - pk00) + t0 .* (pk00 - pk01) + t00 .* (pk01 - pk0)) ./ d;
  b = (t01 .* t01 .* (pk00 - pk0) + t0 .* t0 .* (pk01 - pk00) + t00 .* t00 .* (pk0 - pk01)) ./ d;
  c = (t0 .* t01 .* (t0 - t01) .* pk00 + t01 .* t00 .* (t01 - t00) .* pk0 + t00 .* t0 .* (t00 - t0) .* pk01) ./ d;

tpk =-b./(2.*a);%critical point position

pk =c-b.^2./(4.*a);%critical point value

tpk(isnan(tpk)) = t0(isnan(tpk));
pk(isnan(pk)) = pk0(isnan(pk));
tpk(isnan(pk)) = nan;

A=cat(dim,a,b,c);%parabola coefficients
