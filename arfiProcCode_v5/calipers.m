function D = calipers(ax)
if  nargin == 0
    ax = gca;
end
    [x y] = getline(ax);
    D = sqrt(diff(x)^2+diff(y)^2);