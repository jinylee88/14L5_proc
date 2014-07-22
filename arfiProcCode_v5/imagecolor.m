function [h cb] = imagecolor(cmap,varargin)
ax = gca;
k = get(ax,'Children');
if ~isempty(k)
h0 = imagesc(varargin{:});
set(h0,'visible','off')
if ishold
set(ax,'Children',[k;h0]);
end
end
h = imagesc(varargin{:});
map0 = colormap;
map = colormap(cmap);
N = size(map,1);
colormap(map0)
data = get(h,'CData');
clim = caxis;
idx = min(N,max(1,1+floor((N-1)*((data-clim(1))/diff(clim)))));
data1 = zeros(size(data,1),size(data,2),3);
if nargout == 2
    pos = get(ax,'position');
    cbw = pos(3)*0.05;
    pos1 = pos.*[1 1 0 1]+[0 0 pos(3)-1.5*cbw 0];
    cbpos = pos.*[1 1 0 1] + [pos(3)-cbw 0 cbw 0];
    set(ax,'position',pos1);
    cb = subplot('position',cbpos);
    image(0,linspace(clim(1),clim(2),N),permute(map,[1 3 2]),'parent',cb);
    set(cb,'XTick',[],'YAxisLocation','right','YDir','normal');
    %cbk = get(cb,'Children');
    %set(cbk,'CData',permute(map,[1 3 2]));
end
    
for i = 1:3
    data1(:,:,i) = reshape(map(idx,i),size(data1,1),size(data1,2));
end
set(h,'CData',data1,'CDataMapping','direct');
