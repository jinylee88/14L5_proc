function x = subsampleshift(x0,f0,tau)
if numel(tau) == 1
    x = (ifft(fft(x0)*exp(-1j*2*pi*f0*tau)));
elseif size(tau,1) == 1;
    tau = repmat(tau,size(x0)./size(tau));
    x = (ifft(fft(x0).*exp(-1j*2*pi*f0*tau)));
else
    x = x0;
    for i =1 :size(tau,2)
        for j = 1:size(tau,3);
            for k = 1:size(tau,4);
                x(:,i,j,k) = diag(ifft(fft(x0)*exp(-1j*2*pi*f0*tau(:,i,j,k)).'));
            end
        end
    end
    
end