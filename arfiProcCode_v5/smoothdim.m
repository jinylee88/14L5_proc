function B = smoothdim(A,SPAN,METHOD,DIM)
sz = size(A);
A1 = permute(A,[DIM 1:DIM-1 DIM+1:length(sz)]);
A2 = reshape(A1,size(A1,1),[]);
B2 = A2;
parfor i = 1:size(A2,2);
B2(:,i) = smooth(A2(:,i),SPAN,METHOD);
end
B1 = reshape(B2,size(A1));
B = permute(B1,[2:DIM 1 DIM+1:length(sz)]);
