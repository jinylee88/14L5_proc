function procAll14L5()
resFile = dir('resS1_*.mat');
parFile = dir('parS1_*.mat');
%for i=1:length(resFile)
for i=3:4
timeStamp = resFile(i).name(7:end-4);
prep14L5Data(resFile(i).name, parFile(i).name, [1:7], [9:15]);

[dttps_file] = process_DTTPS(pwd,0);


imData = gen_imData_AdamV1(pwd,5,[10 3],1,7, 1);

save(['imData_' timeStamp '.mat'], 'imData')
end