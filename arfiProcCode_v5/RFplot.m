lidx = 180;
tup = ((0:K-1)*(1/fs))';
tk = tup(lags+lidx);
upmixk = exp(1j*2*pi*f0*tk);
upmix = exp(1j*2*pi*f0*tup);
downmix = exp(-1j*2*pi*f0*tup);

rf1 = real(x1b(:).*upmix(:));
rf2 = real(x2b(:).*upmix(:));
rfk1 = real(x1bt(:,lidx).*upmixk(:));
rfk2 = real(x2bt(:,lidx).*upmixk(:));
clf
subplot(121)
env1 = abs(hilbert(rf1));
env2 = abs(hilbert(rf2));
plot(tup(lidx+lags),env1(lidx+lags),'b-',tup(lidx+lags),env2(lidx+lags),'r-','linewidth',2);
hold all
plot(tup(lidx+lags),-env1(lidx+lags),'b-',tup(lidx+lags),-env2(lidx+lags),'r-','linewidth',2);
plot(tup(lidx+lags),rf1(lidx+lags),'b-',tup(lidx+lags),rf2(lidx+lags),'r-');
%plot(tk,rfk1,'b.--',tk,rfk2,'r.--')
x3b1 = subsampleshift(x2b1(:,lidx),f0,0.25/fs);
x4b1 = subsampleshift(x1b1(:,lidx),f0,-0.25/fs);
rfk3 = real(x3b1.*upmixk(:));
rfk4 = real(x4b1.*upmixk(:));
p3 = plot(tk,rfk3,'b--','linewidth',2);
p4 = plot(tk,rfk4,'r--','linewidth',2);
xlim(tup(lidx)+0.5e-6*[-1 1]);
%axis([0.9e-7 1.15e-7 14 21]);
axis manual
%p3a = plot(tk,abs(hilbert(rfk3)),'b:');
%p4a = plot(tk,abs(hilbert(rfk4)),'b:');
testlags = linspace(-5e-6,5e-6,200)/770;
cc = nan(size(testlags));
subplot(122)
plot(testlags,0*testlags,'k-')
hold on
%plot(-3.5/fs,0,'ks','markersize',16)
pcc = plot(testlags*2,angle(cc),'b-');
pcc2 = plot(testlags*2,abs(cc),'r-');
ylim([-pi/2 pi/2]);
axis manual
plot(testlags*2,0*testlags,'k-')
xlim([min(testlags*2),max(testlags*2)]);
for i = 1:length(testlags);
lagi = sign(testlags(i))*round(abs(testlags(i))*fs);
lagf = testlags(i)-lagi*(1/fs);%sign(testlags(i))*mod(abs(testlags(i)),1/fs);
lidx3 = lidx+lagi;
lidx4 = lidx-lagi;
x3b0 = circshift(x1b,-lagi);
x4b0 = circshift(x2b,lagi);
x3b1 = x3b0(lags+lidx);
x4b1 = x4b0(lags+lidx);
x3b1 = subsampleshift(x3b1,f0,-lagf);
x4b1 = subsampleshift(x4b1,f0,lagf);
tk3 = tup(lags+lidx3);
tk4 = tup(lags+lidx4);
upmixk3 = exp(1j*2*pi*f0*tk3);
upmixk4 = exp(1j*2*pi*f0*tk4);
rfk3 = real(x3b1.*upmixk3(:));
rfk4 = real(x4b1.*upmixk4(:));
downmixk = exp(-1j*2*pi*f0*tk);
x3b2 = hilbert(rfk3).*downmixk;
x4b2 = hilbert(rfk4).*downmixk;
set(p3,'ydata',rfk3);
set(p4,'ydata',rfk4);
%set(p3a,'ydata',abs(hilbert(rfk3)));
%set(p4a,'ydata',abs(hilbert(rfk4)));
cc(i) = x4b2'*x3b2;
set(pcc,'ydata',angle(cc));
set(pcc2,'ydata',exp(-50*(1-abs(cc)/max(abs(cc)))));
pause(0.01);
end