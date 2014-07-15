function [u, Iup, Qup, fdemest] = runLoupas(I, Q, interpFactor, kernelLength, axial, par)
% function [u Iup Qup] = runLoupas(I, Q, interpFactor, axial, par)
%
% Inputs: I - in-phase data
%         Q - quadrature data
%         interpFactor - upsampling factor
%         axial - axial vector (used for demodulation)
%         par - parameters structure generated by arfi_image
%
% Outputs: u - displacement matrix
%          Iup - upsampled in-phase data
%          Qup - upsampled quadrature data

% Setup parameters
fs = par.fs*1e6;
fc = par.fc;
c = par.c; % m/s
kasai_scale = c/(2*pi);

D = size(I); 
D(1) = D(1).*interpFactor;
tstart = tic;
% Iup = zeros(D);
% Qup = zeros(D);
% for i = 1:size(I,2)
%     [tmp1 tmp2] = computeUpsampledIQdata(I(:,i,:),Q(:,i,:),interpFactor);
%     Iup(:,i,:) = reshape(tmp1, D(1), D(3));
%     Qup(:,i,:) = reshape(tmp2, D(1), D(3));
%     for j = 1:size(I,3)
%         [Iup(:,i,j),Qup(:,i,j)] = computeUpsampledIQdata(I(:,i,j),Q(:,i,j),interpFactor);
%     end
% end
[Iup, Qup] = computeUpsampledIQdata(I,Q,interpFactor);
Iup = reshape(Iup, D);
Qup = reshape(Qup, D);
tend = toc(tstart);
fprintf(1, 'Upsampling Computation Time: %0.2fs\n', tend);

fs = fs*interpFactor;
kasai_avg = round(kernelLength*fs/fc);

fdem = par.Apl3.Mod(1).DsF.data*1e6; % frequency dataset values (Hz shift)
frange = par.Apl3.Mod(1).DsF.rr;  % reference ranges (mm)
fc_vec = reshape(interp1(frange, fdem, axial), size(Iup,1), 1);
fdem_vec = fc_vec./fs;

%Compute Displacements
tstart = tic;
u = zeros(size(Iup));
fdemest = zeros(size(Iup));
% N = size(Iup,1);
% parfor i = 1:size(Iup,2)
%     Iref = squeeze(Iup(:,i,1));
%     Qref = squeeze(Qup(:,i,1));
%     u(:,i,:) = loupasParallel(Iref,Qref,squeeze(Iup(:,i,:)),squeeze(Qup(:,i,:)),N,kasai_avg,fdem_vec,fc_vec,kasai_scale);
% end
for i = 1:size(Iup,2)
    if par.ref_idx ~=-1
        Iref = squeeze(Iup(:,i,par.ref_idx));
        Qref = squeeze(Qup(:,i,par.ref_idx));
    end
%     u(:,i,:) = loupasParallel(Iref,Qref,squeeze(Iup(:,i,:)),squeeze(Qup(:,i,:)),N,kasai_avg,fdem_vec,fc_vec,kasai_scale);
    if i==1
        fprintf(1, 'Displacement Estimation for Beam %d/%d', i, size(Iup,2));
    elseif i==size(Iup,2)
        tmpS = sprintf('%d/%d', i-1, size(Iup,2));
        fprintf(1, repmat('\b', [1 length(tmpS)]));
        fprintf(1, '%d/%d', i, size(Iup,2));
        fprintf(1, '\n');
    else
        tmpS = sprintf('%d/%d', i-1, size(Iup,2));
        fprintf(1, repmat('\b', [1 length(tmpS)]));
        fprintf(1, '%d/%d', i, size(Iup,2));
    end
    if par.ref_idx ==-1
        for k = 1:size(Iup,3)-1
            Iref = squeeze(Iup(:,i,k));
            Qref = squeeze(Qup(:,i,k));
            Idisp = squeeze(Iup(:,i,k+1));
            Qdisp = squeeze(Qup(:,i,k+1));
            u_inc(:,i,k) = loupas(Iref,Qref,Idisp,Qdisp,size(Iup,1),kasai_avg,fdem_vec,fc_vec,kasai_scale);
        end
    else
        for k = 1:size(Iup,3)
            Idisp = squeeze(Iup(:,i,k));
            Qdisp = squeeze(Qup(:,i,k));
            u(:,i,k) = loupas(Iref,Qref,Idisp,Qdisp,size(Iup,1),kasai_avg,fdem_vec,fc_vec,kasai_scale);
            %         [u(:,i,k), fdemest(:,i,k)] = loupas(Iref,Qref,Idisp,Qdisp,size(Iup,1),kasai_avg,fdem_vec,fc_vec,kasai_scale,out);
        end
    end
end

% Generating accumulated displacements from incremental
if par.ref_idx ==-1
    u(:,:,par.nref-1:-1:1) = cumsum(u_inc(:,:,par.nref-1:-1:1),3);
    u(:,:,par.nref+(1:par.npush+2)) = u_inc(:,:,par.nref+(0:par.npush+1));
    u(:,:,par.nref+par.npush+3:end) = cumsum(u_inc(:,:,par.nref+par.npush+2:end),3);
end
    u = -u.*1e6;
    u = u(1:end-kasai_avg-1, :, :, :);
tend = toc(tstart);

fprintf(1, 'Displacement Computation Time: %0.2fs\n', tend);