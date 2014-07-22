  function [v] = computeFluxSamTrack(data,t,k_length,c,f0)

 [L,W,H] = size(data);
 data0 = reshape(data(:,:,1:end-1),L,W*(H-1));
 data1 = reshape(data(:,:,2:end),L,W*(H-1));
 dt01 = diff(t);
 dt10 = diff(t(end:-1:1));
 swin = ones(k_length,1);  
 v = kasai(data0,data1,dt01,swin,c,f0,L,W,H);
% v10 = kasai(data1,data0,-dt01,swin,c,f0,L,W,H);
  
  
function v = kasai(data0,data1,dt,swin,c,f0,L,W,H)
 nu01 = imag(data1).*real(data0)-real(data1).*imag(data0);
 de01 = real(data1).*real(data0)+imag(data1).*imag(data0);
 nu01 = conv2(nu01,swin,'same');
 de01 = conv2(de01,swin,'same');
 u  = atan2(nu01,de01);
 u  = u * (c/(f0*4*pi)); 
u = reshape(u,L,W,H-1);
dt = repmat(permute(dt(:),[2 3 1]),L,W);
v = u./dt;
return

 