function procArfi_wrapper(datadir,i);
addpath /getlab/pjh7/SC2000/arfiProcCode_v4/
cd(datadir);
fNames = dir('SWIF_AData*.bin');
procArfi_v4(fNames(2*(i-1)+1).name,'method','PesaventoFlux','kernellength',9,'removeKickBack',0,'interpFactor',1);
