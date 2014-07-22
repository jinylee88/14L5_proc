
function TTPdynWrapper(fileidx)
    warning off all
    filepaths = {'/getlab/pjh7/SC2000/data/20130530_LV_Mmode/','/getlab/pjh7/SC2000/data/20130530_RV_Mmode/','/getlab/pjh7/SC2000/data/20130530_RA_Mmode/'};
    tic
    D = {};
    kz = 5;
    kt = 2;
    bandpass = [100 1000];
    for filepathidx = 1:3
    filepath = filepaths{filepathidx};
    method = 'PesaventoFlux';
    respath = fullfile(filepath,'res',method);
    d = dir(fullfile(respath,'res_all_*.mat'));
    for i = 1:length(d);
    D{end+1} = fullfile(respath,d(i).name);
    end
    end
    resfile = D{fileidx};
    timestamp = RetrieveTimeStamp(resfile);
    fprintf('loading %s...',resfile);
    res = load(resfile);
    fprintf('done\n');
    fprintf('differentiating...')
    v0 = diff(double(res.arfidata0_all)*res.arfi_scale,1,3);
    v0(:,:,56,1:end-1) = 0.5*(v0(:,:,55,1:end-1)+v0(:,:,1,2:end));
    v0 = v0(:,:,:);
    v0 = v0(:,:,[end 1:end-1]);
    [v2 t1] = differentiateDisplacements(cumsum(v0,3),res.t,bandpass);
    v2 = reshape(v2(:,:,[1:end end]),size(res.arfidata0_all));
    v2 = v2(:,:,1:end-1,:);
    clear v0
    fprintf('done\n')
    fprintf('Left-Right filtering...')
    [vL vR] = LeftRightFilter(v2-repmat(median(v2,2),[1 32 1 1]));
    clear v2
    vLR = [vL(:,1:16,:,:) vR(:,17:32,:,:)];
    clear vL vR
    fprintf('done\n');
    [DX DZ] = sectorCoordinates(res.theta,res.theta,res.axial,res.apex);
 
    for i = 1:size(vLR,1);
        fprintf('Axial location %0.0f/%0.0f...',i,size(vLR,1));
        zidx = max(1,i-kz):min(size(vLR,1),i+kz);
        for j = 1:192;
            tidx = max(1,j-kt):min(size(vLR,4),j+kt);
            dxdt(i,j) = TTPdyn(t1,mean(DX(zidx,:),1),squeeze(median(mean(vLR(zidx,:,:,tidx),1),4)));
        end;
    end
    outdir = fullfile(filepath,'dxdt','dynttp');
    if ~exist(outdir,'dir')
        mkdir(outdir)
    end
    outfile = fullfile(outdir,sprintf('dxdt_%s.mat',timestamp));
    save(outfile,'dxdt','kz','kt','bandpass')