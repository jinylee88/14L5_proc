  function [arfidata] = computeDisplacementsSamTrack(data,preTracks,k_length,c,f0)
%Run Kasai's algorithm
    kdata = data;%permute(data,[1 3 2]);
    kdata = kdata(:,:,preTracks+1:end); 
    [L,W,H] = size(kdata);
    swin=ones(k_length,1)/k_length;
    dref = repmat(kdata(:,1:W),1,H-1);
    kdata = kdata(:,W+1:end);
    nu=imag(kdata).*real(dref)-real(kdata).*imag(dref);
    de=real(kdata).*real(dref)+imag(kdata).*imag(dref);
    nu=conv2(nu,swin,'same');
    de=conv2(de,swin,'same');
    displace_est=atan2(nu,de);
    displace_est = displace_est * (c/(f0*4*pi)); %um
    arfidata_kasai = reshape(displace_est,L,W,H-1);
    %Track pretracks in reverse
    if preTracks>0;
    kdata1 = data;%permute(data,[1 3 2]);
    kdata1 = kdata1(:,:,[preTracks+1 preTracks+1:-1:1]);
    [L,W,H] = size(kdata1);
    dref = repmat(kdata1(:,1:W),1,H-1);
    kdata1 = kdata1(:,W+1:end);
    nu=imag(kdata1).*real(dref)-real(kdata1).*imag(dref);
    de=real(kdata1).*real(dref)+imag(kdata1).*imag(dref);
    nu=conv2(nu,swin,'same');
    de=conv2(de,swin,'same');
    displace_est=atan2(nu,de);
    displace_est = displace_est * (c/(f0*4*pi)); %um
    pretrack = reshape(displace_est,L,W,H-1);
    arfidata_kasai = cat(3,pretrack(:,:,end:-1:1),arfidata_kasai);
    else
    arfidata_kasai = cat(3,zeros(size(arfidata_kasai,1),size(arfidata_kasai,2),1),arfidata_kasai);
    end
    arfidata = arfidata_kasai;
    