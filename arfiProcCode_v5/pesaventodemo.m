x1b0 = IQ0(:,8);
Tau = 0.02e-6+(rand(size(x1b))-0.5)*1e-8;
Tau = filter2(ones(7,1)/7,Tau);
Tau = 1e-8;
x1b = subsampleshift(x1b0,f0,-Tau/2);
x2b = subsampleshift(x1b0,f0,Tau/2);
%x2b = IQ0(:,11);
N = 20;
Tk0 = 0;
Tk = 0;
f0 = 6e6;
w0 = 2*pi*f0;
wm = w0;
t = (0:length(x1b)-1)*(1/(par.fs*1e6));
Ts = mean(diff(t));
fs = 1/Ts;
Tw = 2e-4/770;
dt = 0.1e-6/770;
lags = (-Tw/2):dt:(Tw/2);
phz = conj(x1b).*(x2b);
K = ones(ceil(Tw/Ts),1);
K = ones(1,1,1);
I = Ts*ones(size(x1b));
Denom = convn(I,K,'same');
X2bt = (subsampleshift(x2b,f0,-lags/2));
X1btc = conj(subsampleshift(x1b,f0,lags/2));
tup = 0:(1/80e6):t(end);
x1bup = interp1(t,x1b,tup,'spline');
x1rf =  x1bup.*exp(1j*w0*tup);
x2bup = interp1(t,x2b,tup,'spline');
x2rf = x2bup.*exp(1j*w0*tup);

tau0 = zeros(length(x1b),1);
tau = tau0;   
for n = 1:N
    for ii = linspace(0,1,50)
    x1bt = subsampleshift(x1b,f0,-(ii*tau + (1-ii)*tau0)/2);
    x2bt = subsampleshift(x2b,f0,(ii*tau + (1-ii)*tau0)/2);
    x1rft = (interp1(t,x1bt,tup,'spline').*exp(1j*w0*tup));
    x2rft = (interp1(t,x2bt,tup,'spline').*exp(1j*w0*tup));
    plot(tup,real(x1rf),'b.:',tup,real(x2rf),'r.:',...
        tup,real(x1rft),'k-',tup,real(x2rft),'m-',...
        tup,abs(hilbert(real(x1rf))),'b-',...
        tup,abs(hilbert(real(x1rft))),'k-',...
        tup,-1*abs(hilbert(real(x2rf))),'r-',...
        tup,-1*abs(hilbert(real(x2rft))),'m-'...
        );xlim([6 11]*1e-6);ylim([-250 250]);
    pause(0.01);
    drawnow;
    end
    tau0 = tau;
    corrcoeff = (convn(conj(x1bt).*x2bt,K,'same')).*Denom;
    phzcorrect = angle(exp(-1j*wm*tau0).*(corrcoeff));
    tau = tau0 + (1/(w0))*phzcorrect;    
    tau = -1*tau;
    %set(p,'xdata',tau);
    %subplot(222)
   
%     plot(real(x1rf),tup,real(x2rf),tup);
%     axis ij
%     ylim([5e-6 6e-6]);
%     subplot(224)
%     phz0 = angle(conj(x1b).*x2b);
%     newphz = angle(conj(x1b).*subsampleshift(x2b,f0,-tau));
%     set(p3,'xdata',tau(20),'ydata',newphz(20));
    pause;
end




