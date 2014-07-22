function dxdt = TTPdyn(t,dx,slice)
disp_on = 0;
dx = round(1e6*dx(:))/1e6;
t = t(:);
%centeridx = find(abs(dx)==min(abs(dx)));
dx0 = 1;
firstkernel = 1;
kernelstep = 0.5;
maxdx = 4;
slicetmp = slice;
if disp_on
clf
im = imagesc(t,dx,slicetmp,[-10 10]);
hold on
p = plot(0,0,'wo-');
p1 = plot(0,0,'w.');
end
%for i = 1:length(dx);
%    [pks{i} locs{i}] = findpeaks(slice(i,:),'minpeakdistance',5,'minpeakheight',max(slice(i,:))/4);
%    dxs{i} = repmat(dx(i),size(pks{i}));
%end
%     F = fit(abs([dxs{:}]'),t([locs{:}]'),'poly1','weights',[pks{:}],'robust','on');
%     dxdt = 1./F.p1;
%     if disp_on
%         set(p,'xdata',F.p1*abs(dx)+F.p2,'ydata',dx);
%         set(p1,'xdata',t([locs{:}]),'ydata',[dxs{:}]);
%     end
[T DX] = meshgrid(t,dx);
for dx1 = dx0+firstkernel:kernelstep:maxdx
    xidx = (abs(dx)>=dx0) & (abs(dx)<=dx1);
    %F = fit(abs([dxs{xidx}]'),t([locs{xidx}]'),'poly1','weights',[pks{xidx}],'robust','on')
    %P = [F.p2 F.p1];
    [pk tpk] = subsamplepeak(t,slicetmp,2);
    P = polyfit(abs(dx(xidx)),tpk(xidx),1);
    tfit = polyval(P,abs(dx));
    if disp_on
    set(im,'cdata',slicetmp);
    set(p,'xdata',tfit(xidx),'ydata',dx(xidx));
    set(p1,'xdata',tpk(xidx),'ydata',dx(xidx));
    pause(0.1);
    end
    if 1/P(1)<20 && 1/P(1)>0 && P(2)>-1
    slicetmp = slice;
    slicetmp(((abs(DX)./(T-P(2)+0.5))>(1/P(1)+1))&((T-P(2)+0.5)>=0)) = nan;
    slicetmp(((abs(DX)./(T-P(2)-1))<(1/P(1)-0.5))&((T-P(2)-1)>=0)) = nan;
    %for i = 1:size(slicetmp,1)
    %    slicetmp(i,:) = slice(i,:);
    %    win = max(1,(abs(dx(i))-1)/3);
    %    slicetmp(i,t<tfit(i)-win | t>tfit(i)+1.5*win) = nan;
    %end
    end
end
dxdt = 1/P(1);    
