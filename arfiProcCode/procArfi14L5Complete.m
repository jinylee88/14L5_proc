function procArfi14L5Complete(fname, parFile, interpFactor, kernelLength, ccmode, ref_idx)

%This code will process a file acquired using alternating focused/unfocused
%tracking beams for the combined SWEI/ARFI acquisitions. The output will be
%a resA_* file which is the arfi data and a resS_* file which is the SWEI
%data, although the SWEI data at this stage will still have all track
%locations including the on axis ones.
if nargin==0
    procArfi14L5ModFocusedARFI();
    procArfi14L5ModUnfocusedSWEI();
elseif nargin==1
    procArfi14L5ModFocusedARFI(fname);
    procArfi14L5ModUnfocusedSWEI(fname);
elseif nargin==2
    procArfi14L5ModFocusedARFI(fname, parFile);
    procArfi14L5ModUnfocusedSWEI(fname, parFile);
elseif nargin==3
    procArfi14L5ModFocusedARFI(fname, parFile,interpFactor);
    procArfi14L5ModUnfocusedSWEI(fname, parFile,interpFactor);
elseif nargin==4
    procArfi14L5ModFocusedARFI(fname, parFile,interpFactor,kernelLength);
    procArfi14L5ModUnfocusedSWEI(fname, parFile,interpFactor,kernelLength);
elseif nargin==5
    procArfi14L5ModFocusedARFI(fname, parFile,interpFactor,kernelLength, ccmode);
    procArfi14L5ModUnfocusedSWEI(fname, parFile,interpFactor,kernelLength, ccmode);
elseif nargin==6
    procArfi14L5ModFocusedARFI(fname, parFile,interpFactor,kernelLength, ccmode, ref_idx);
    procArfi14L5ModUnfocusedSWEI(fname, parFile,interpFactor,kernelLength, ccmode, ref_idx);
else
    procArfi14L5ModFocusedARFI();
    procArfi14L5ModUnfocusedSWEI();
end
end