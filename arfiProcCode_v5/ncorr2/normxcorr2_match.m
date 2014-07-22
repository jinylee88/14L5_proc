function C = normxcorr2_match(template,img)
img_padded=padarray(img,[size(template,1)-1,size(template,2)-1]);
C=normxcorr2_mex(template, img_padded, 'valid');