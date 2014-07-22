function H = surfStack(Z,X,Y,V,A,dim,DownSample);
%SURFSTACK show volume of data as stack of surfaces
% H = surfStack(Z,X,Y,V,A,dim,DownSample)
% surfStack takes the volumetric data specified in V and the transparency
% (alpha) data stored in A, in displays it in 3-D according to the axis
% vectors Z,X,and Y, corresponding to the 1st, 2nd, and 3rd dimensions of
% V.
%
% INPUTS:
%   Z: axis vector in the first (axial) direction (Mx1)
%   X: axis vector in the second (lateral) direction (Nx1)
%   Y: axis vector in the third (elevational) direction (Px1)
%   V: volumetric data, scaled (MxNxP) or RGB (MxNxPx3)
%   A: transparency data, scaled {0,1}, (MxNxP).
%   dim: data dimension normal to each plane (1,2, or 3) - defaults to 3
%   DownSample: downsampling factor applying to entire volume 
%       Three element vector [ds_z ds_z ds_y] specifies per dimension, and
%       a single value will specifiy all three dimensions
%
% OUTPUTS:
%   H: Handle vector pointing to the surfaces
%
% BE AWARE THAT RENDERING TIMES CAN BE VERY LONG FOR LARGE DATASETS!
% USING SEMI-TRANSPARENCY (0<VALUE<1) ALSO INCREASES RENDERING TIME!
%
% Use 'cameratoolbar show' to get more advanced camera control to specify
% a new CameraUpVector for easier rotation of complex volumes
%
% See also CAMERATOOLBAR

% Version History
% Created 09/02/2013
% Peter Hollender

if ~exist('dim','var')
    dim = 3;
end
if ~exist('DownSample','var')
    DownSample = [1 1 1];
end
if length(DownSample) == 1
    DownSample = DownSample*[1 1 1];
end
holdState = get(gca,'nextPlot');
set(gca,'nextPlot','add');
[sz(1) sz(2) sz(3) rgbn] = size(V);
[sz1(1) sz1(2) sz1(3)] = size(A);
A = repmat(A,[sz./sz1]);
z = Z(1:DownSample(1):end);
y = Y(1:DownSample(3):end);
x = X(1:DownSample(2):end);
v = V(1:DownSample(1):end,1:DownSample(2):end,1:DownSample(3):end,:);
a = A(1:DownSample(1):end,1:DownSample(2):end,1:DownSample(3):end);
dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));
switch dim
        case 1
        xdata = [x(1)-2*dx;x(1)-dx;x(:);x(end)+dx] - dx/2;
        ydata = [y(1)-2*dy;y(1)-dy;y(:);y(end)+dy] - dy/2;
        [xdata ydata] = meshgrid(xdata,ydata);
        zdata = z;
        for i = 1:length(z)
        cdata =  padarray(padarray(...
            permute(squeeze(v(i,:,:,:)),[2 1 3]),...
            [2 2 0],'pre'),[1 1 0],'post');
        adata = padarray(padarray(...
            permute(squeeze(double(a(i,:,:))),[2 1 3]),...
            [2 2 0],'post'),[1 1 0],'pre');
        H(i) = surf(xdata,ydata,zdata(i)*ones(size(ydata)),cdata,...
            'edgecolor','none','alphadata',adata,...
            'facealpha','flat','alphadatamapping','none');
        end
        case 2
        zdata = [z(1)-2*dz;z(1)-dz;z(:);z(end)+dz] - dz/2;
        ydata = [y(1)-2*dy;y(1)-dy;y(:);y(end)+dy] - dy/2;
        [zdata ydata] = meshgrid(zdata,ydata);
        xdata = x;
        for i = 1:length(x)
        cdata =  padarray(padarray(...
            permute(squeeze(v(:,i,:,:)),[2 1 3]),...
            [2 2 0],'pre'),[1 1 0],'post');
        adata = padarray(padarray(...
            permute(squeeze(double(a(:,i,:))),[2 1 3]),...
            [2 2 0],'post'),[1 1 0],'pre');
        H(i) = surf(xdata(i)*ones(size(ydata)),ydata,zdata,cdata,...
            'edgecolor','none','alphadata',adata,...
            'facealpha','flat','alphadatamapping','none');
        end   
        
    case 3
        xdata = [x(1)-2*dx;x(1)-dx;x(:);x(end)+dx] - dx/2;
        zdata = [z(1)-2*dz;z(1)-dz;z(:);z(end)+dz] - dz/2;
        [xdata zdata] = meshgrid(xdata,zdata);
        ydata = y;
        for i = 1:length(y)
        cdata =  padarray(padarray(...
            permute(squeeze(v(:,:,i,:)),[1 2 3]),...
            [2 2 0],'pre'),[1 1 0],'post');
        adata = padarray(padarray(...
            permute(squeeze(double(a(:,:,i))),[1 2 3]),...
            [2 2 0],'post'),[1 1 0],'pre');
        H(i) = surf(xdata,ydata(i)*ones(size(zdata)),zdata,cdata,...
            'edgecolor','none','alphadata',adata,...
            'facealpha','flat','alphadatamapping','none');
        end

end
set(gca,'nextPlot',holdState);