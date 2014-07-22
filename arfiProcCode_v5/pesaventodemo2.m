function Tau = Pesavento(x1b,x2b,Klen,N)
if ~exist('Klen','var')
Klen = 1;
end
if ~exist('N','var')
N = 3;
end
tau = 0;
for k = 1:size(x1b,1);
    kidx = max(1,k-floor(Klen/2)):min(size(x1b,1),k+floor(Klen/2));
    tup1 = tup(tup>=min(t(kidx)) & tup<=max(t(kidx)));
    %x1bup = interp1(t(kidx),x1b(kidx),tup1,'spline');
    %x2bup = interp1(t(kidx),x2b(kidx),tup1,'spline');
    %x1rf =  x1bup.*exp(1j*w0*tup1);
    %x2rf = x2bup.*exp(1j*w0*tup1);
    for n = 1:N
    x1bt = subsampleshift(x1b(kidx),f0,tau/2);
    x2bt = subsampleshift(x2b(kidx),f0,-tau/2);
    %x1tup = interp1(t(kidx),x1bt,tup1,'spline');
    %x2tup = interp1(t(kidx),x2bt,tup1,'spline');
    %x1rft =  x1tup.*exp(1j*w0*tup1);
    %x2rft =  x2tup.*exp(1j*w0*tup1);
    corrVal = x2bt(:)'*x1bt(:);
    %corrValrf = x2rft(:)'*x1rft(:);
    corrValPhasCrct = angle(exp(-1j*wm*(tau)/fs)*corrVal);
    tau = tau + corrValPhasCrct/wm;
    end
    Tau(k) = tau;
end
    