
function CT = sws_SMURF(TTP,DX,n,par)
            [ii jj] = sort(par.pushbeamNum);
            B = [-1*ones(1,ceil(n/2)) zeros(1,mod(n+1,2)) ones(1,ceil(n/2))];
            B = permute(B,[1 3 2]);
            CT = squeeze(nanmedian(convn(DX(:,:,jj),B,'same')./convn(TTP(:,:,jj).*sign(DX(:,:,jj)),B,'same'),2));
end
