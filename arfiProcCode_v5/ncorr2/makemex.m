% Run this code to build the mex function.
cmd = 'mex ';
cputype = computer('arch');
if strcmp(cputype(end-1:end),'64')
    cmd = [cmd '-v -largeArrayDims -lrt '];
end
d = dir('cv_src/*.cpp');
cmd = [cmd 'normxcorr2_mex.cpp' sprintf(' cv_src/%s',d(:).name)];
disp(cmd)
eval(cmd)

clear cmd cputype