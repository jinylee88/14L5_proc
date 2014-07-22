function [Lags x2] = PesaventoGross(x1b,x2b,Klen,SrchLen);
Lag = -SrchLen:SrchLen;
x2 = x2b;
K = size(x1b,1);
  for k = 1:K;
    kidx = max(1,k-floor(Klen/2)):min(K,k+floor(Klen/2));
    srchidx = min(K,max(1,(kidx(1)-SrchLen):(kidx(end)+SrchLen)));
    grossCorr = ncorrAccel_V3(abs(x1b(kidx,:,:)),abs(x2b(srchidx,:,:)));
    [pkCorr,Idx(k,:,:)] = max(grossCorr);
  end
  for j = 1:size(Idx,3);
      Idx(:,:,j) = round(medfilt2(Idx(:,:,j),[Klen,1]));
  end
  Idx1 = repmat(permute(1:K,[2 1 3]),[1 size(x1b,2) size(x1b,3)]);
  Idx2 = max(1,min(K,Idx1+Lag(Idx)));
  Lags = Idx2-Idx1;
    for i = 1:size(x1b,2);
      for j = 1:size(x1b,3)
          x2(:,i,j) = x2b(Idx2(:,i,j),i,j);
      end
    end
 