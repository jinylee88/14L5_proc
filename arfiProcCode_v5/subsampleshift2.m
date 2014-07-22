function x = subsampleshift2(x0,f0,fs,shift)
shifti = sign(shift).*round(abs(shift)*fs);
shiftfi = shift*fs - shifti;
shiftt = shiftfi/fs;
if numel(shift) == 1
    x = (ifft(fft(circshift(x1,shifti))*exp(-1j*2*pi*f0*shiftt)));
else
    shift = repmat(shift,[size(x0,1)./size(shift,1) 1 1 1]);
    x = x0;
    for i =1 :size(shift,2)
        for j = 1:size(shift,3);
            for k = 1:size(shift,4);
                x(:,i,j,k) = diag(ifft(fft(circshift(x0(:,i,j,k),shifti())*exp(-1j*2*pi*f0*shiftt(:,i,j,k)).'));
            end
        end
    end
    
end