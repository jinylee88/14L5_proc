function F = dft(x,y,f)
F = zeros(size(f));
for k = 1:length(f)
F(k) = sum(y.*exp(-1j*2*pi*f(k)*x));
end
%F = exp(-1j*2*pi*(f(:)*x(:)'))*y(:)