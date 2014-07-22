
figure(figHandles.IQBMode);clf
set(figHandles.IQBMode,'Name','IQ B-Mode','NumberTitle','off')
image(blat,baxial,bimg(:,:,[1 1 1]));
xlabel('x (mm)')
ylabel('z (mm)')
title(timeStamp)
colormap gray
colorbar('location','northoutside')
axis image
axis(rt.axisLimits);

switch swifParams.imagingModeShort
    case 'mmode'
        arfidata = squeeze(arfidata);
        figure(figHandles.rtSlices);
        if frmii == 1
            clf;
        end
        rowmax = min(16,length(Frameindices));
        ii = mod(frmii-1,rowmax)+1;
        jj = ceil(frmii/rowmax);
        colmax = ceil(length(Frameindices)/rowmax);
        sp = subplot('position',[(ii-1)*(1/rowmax) 1-(jj*(1/colmax)) (1/rowmax) (1/colmax)]);
        axidx = find(axial>mean(par.pushFocalDepth)-2 & axial<mean(par.pushFocalDepth));
        slice = diff(squeeze(mean(arfidata(axidx,:,:))),1,2);
        slice = slice-repmat(mean(slice),size(slice,1),1);
        imagesc(1:size(arfidata,3),lat(1,:),slice,'Parent',sp);caxis(sp,[-1 1]);
        set(sp,'YTick',[],'XTick',[]);
        colormap(Hollender(64));
        
        if flags.displayflag > 1
        figure(figHandles.rtMov);clf
        A = par.trackParams.specList.scanSpecs.apexMm(3);
        [TH R] = meshgrid(double(lat(1,:)),double(axial));
        X = R.*sind(TH) - A*tand(TH);
        Y = R.*cosd(TH);
        im = surf(X,Y,0*R,double(arfidata(:,:,1)));view([0 90]);set(im,'edgecolor','none');axis ij;axis image;caxis([0 5]);ylim([0 30]);colormap jet
        hold on
        bim2 = image(blat,baxial,bimg(:,:,[1 1 1]));
        set(bim2,'alphadata',bimg(:,:,1)<3/255);
        for i = 1:size(arfidata,3);set(im,'cdata',double(arfidata(:,:,i)));tic;drawnow;while(toc<0.05);end;end
        end
        
    case 'swei'
        figure(figHandles.rtMov);clf
        set(figHandles.rtMov,'Name','Real-Time','NumberTitle','off')
        im = imagesc(lat,axial,squeeze(arfidata(:,:,1)),rt.dispLim);
        set(im,'Alphadata',2*min(0,cc(:,:,1) - 0.5));
        set(gca,'color','k')
        xlabel('lat')
        ylabel('axial (mm)');
        for i = 2:length(t);
            set(im,'CData',squeeze(arfidata(:,:,i)),'Alphadata',5*max(0,double(squeeze(cc(:,:,i)))/255 - 0.8));
            pause(0.01);
        end
        [aimg aaxial alat] = Scan_Convert_Sector(squeeze(arfidata),scB);
        [cc_onax] = Scan_Convert_Sector(squeeze(cc),scB);
        
    case 'arfiswei'
        figure(figHandles.rtMov);clf
        set(figHandles.rtMov,'Name','Real-Time','NumberTitle','off')
        hold on
        set(gca,'color','k');
        axis([min(lat(:)) max(lat(:)) min(axial) max(axial)])
        axis ij
        
        BkH = image([min(lat(1,:))-1 max(lat(1,:))+1],axial,zeros(size(arfidata,1),size(arfidata,2),3));
        xlabel('angle (deg)')
        xlabel('z (mm)');
        cb = colorbar('location','northoutside');
        ylabel(cb,'displacement (\mum)')
        M = [];
        if flags.displayflag>1
            for j = 1:size(arfidata,3)
                set(BkH,'XData',[min(lat(j,:))-1 max(lat(j,:))+1])
                k = get(gca,'Children');
                set(gca,'Children',[BkH;k(k~=BkH)])
                Rh = imagesc(lat(j,Ridx),axial,squeeze(arfidata(:,Ridx,j,1)),[rt.dispLim]);
                set(Rh,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Ridx,j,1)) - rt.ccRng(1)));
                Ch = imagesc(lat(j,Cidx),axial,squeeze(arfidata(:,Cidx,j,1)),[rt.dispLim]);
                set(Ch,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Cidx,j,1)) - rt.ccRng(1)));
                Lh = imagesc(lat(j,Lidx),axial,squeeze(arfidata(:,Lidx,j,1)),[rt.dispLim]);
                set(Lh,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Lidx,j,1)) - rt.ccRng(1)));
                drawnow
                for i = 2:size(arfidata,4);
                    tic
                    set(Rh,'CData',squeeze(arfidata(:,Ridx,j,i)))
                    set(Rh,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Ridx,j,i)) - rt.ccRng(1)));
                    set(Ch,'CData',squeeze(arfidata(:,Cidx,j,i)));
                    set(Ch,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Cidx,j,i)) - rt.ccRng(1)));
                    set(Lh,'CData',squeeze(arfidata(:,Lidx,j,i)));
                    set(Lh,'Alphadata',(1/diff(rt.ccRng))*max(0,squeeze(cc(:,Lidx,j,i)) - rt.ccRng(1)));
                    drawnow
                    while(toc<0.5*diff(t(i+[-1 0])));
                    end
                end
                delete(Lh);
                delete(Rh);
                set(Ch,'CData',squeeze(arfidata(:,Cidx,j,par.nref+4)))
                set(Ch,'Alphadata',(1/diff(rt.ccRng))*max(0,mean(squeeze(cc(:,Cidx,j,par.nref+[8:12])),3) - rt.ccRng(1)));
                drawnow
            end
            delete(BkH);
        end
        aonax = reshape(arfidata(:,Cidx,:,:),[size(arfidata,1),length(Cidx)*par.numBeamGroups,length(t)]);
        conax = reshape(cc(:,Cidx,:,:),[size(arfidata,1),length(Cidx)*par.numBeamGroups,length(t)]);
        thonax = reshape(lat(:,Cidx),1,[]);
        sconax = struct(...
            'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
            'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
            'latinc',.1e-3,...
            'axialmin',-1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'axialmax',1e-2*(par.imagingdepth)- 1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'axialinc',.1e-3,...
            'min_phi',min(thonax),...
            'span_phi',max(thonax)-min(thonax),...
            'apex',0.1*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'fsiq',par.fs*dispEst.dispEst.interpFactor*1e6...
            );
        [arfidata_onax axial_onax lat_onax ] = Scan_Convert_Sector(aonax,sconax);
        [cc_onax] = Scan_Convert_Sector(conax,sconax);
        
        figure(figHandles.ARFI);clf
        im = imagesc(lat_onax,axial_onax,arfidata_onax(:,:,par.nref+4),rt.dispLim);
        hold on
        im0 = imagesc(lat_onax,axial_onax,arfidata_onax(:,:,par.nref+4),rt.dispLim);
        set(im,'Alphadata',(1/diff(rt.ccRng))*max(0,mean(cc_onax(:,:,1:par.nref),3) - rt.ccRng(1)));
        set(gca,'color','k')
        xlabel('x (mm)')
        ylabel('z (mm)')
        cb = colorbar('northoutside');
        xlabel(cb,'displacement \mum')
        bim = image(blat,baxial,bimg(:,:,[1 1 1]));
        set(gca,'Children',[im bim im0])
        axis image
        axis(rt.axisLimits)
        set(figHandles.ARFI,'Name','ARFI','NumberTitle','off');
        title(sprintf('Frame #%03.0f (%s)',frmidx,timeStamp))
        print('-dpng',fullfile(outputdatadir,sprintf('%s_ARFI_%03.0f',timeStamp,frmidx)));
        
    case 'arfi'

        if strcmpi(par.mode,'mmode')
            ii = reshape(reshape(1:par.numBeamGroups,[],2)',[],1);
            ARFIDATA1 = reshape(arfidata,length(axial),par.numBeamGroups,par.nBeams,length(t));
            for i = 1:par.numBeamGroups;
                arfidata1(:,ii(i),:) = mean(ARFIDATA1(:,ii(i),i,:),2);end
            clear ARFIDATA1;
            [TH R] = meshgrid(double(unique(lat)),double(axial));
        else
            arfidata1 = arfidata;
            [TH R] = meshgrid(double(lat(1,:)),double(axial));
        end
        A = par.trackParams.specList.scanSpecs.apexMm(3);
        X = R.*sind(TH) - A*tand(TH);
        Y = R.*cosd(TH);
        sconax = struct(...
            'latmin',1e-2*sind(-45)*par.imagingdepth,...  1e-2*sin(min(btheta))*par.imagingdepth,...
            'latmax',1e-2.*sind(45)*par.imagingdepth,... 1e-2*sin(max(btheta))*par.imagingdepth,...
            'latinc',.1e-3,...
            'axialmin',-1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'axialmax',1e-2*(par.imagingdepth)- 1e-3*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'axialinc',.1e-3,...
            'min_phi',min(lat),...
            'span_phi',max(lat)-min(lat),...
            'apex',0.1*par.trackParams.paramList.xdcrParams.FovBMode.VectorApexAzimZMm,...
            'fsiq',par.fs*dispEst.interpFactor*1e6...
            );
        [arfidata_onax axial_onax lat_onax ] = Scan_Convert_Sector(arfidata1,sconax);
        [cc_onax] = Scan_Convert_Sector(cc,sconax);
        figure(figHandles.ARFI);clf
        im0 = imagesc(lat_onax,axial_onax,arfidata_onax(:,:,par.nref+4),rt.dispLim);
        hold on
        if rt.depthCorrect
        %arfidata_fix = arfidata1./repmat(nanmedian(arfidata1,2),[1 size(arfidata1,2) 1]);
        %im = imsurf(unique(lat),axial,A,double(arfidata_fix(:,:,12)));
        R0 = sqrt((Y-7).^2+X.^2);
        normalization = max(1,(R0/2).^0.35);
        %normalization(R0>10) = 0;
        im  = imsurf(unique(lat),axial,A,double(arfidata1(:,:,12).*normalization));
        else
        im = imsurf(unique(lat),axial,A,double(arfidata1(:,:,12)));
        end
        %im = surf(X,Y,0*R,double(max(arfidata1(:,:,9+[0:10]),[],3)));view([0 90]);set(im,'edgecolor','none');axis ij;axis image;
        set(gca,'color','k')
        xlabel('x (mm)')
        ylabel('z (mm)')
        cb = colorbar('NorthOutside');
        xlabel(cb,'displacement (\mum)')
        bim = image(blat,baxial,bimg(:,:,[1 1 1]));
        bim2 = image(blat,baxial,bimg(:,:,[1 1 1]));
        set(bim2,'alphadata',min(0.8,1-bimg(:,:,1).^0.5));
        set(bim2,'alphadata',min(0.9,bimg(:,:,1)<20/255))

        set(gca,'Children',[bim2 im bim im0])
        if ~rt.bMask;set(bim2,'visible','off');end
        axis image
        axis(rt.axisLimits)
        set(figHandles.ARFI,'Name','ARFI','NumberTitle','off');
        drawnow
        
            if strcmpi(par.mode,'mmode')
            arfidata2 = arfidata(:,:,[1:5 8:end]);
            pushbeamNum = reshape(repmat(par.pushbeamNum,par.nBeams,1),1,[]);
            th1 = interp1(1:par.nBeams,lat(1:par.nBeams),pushbeamNum);
            dth = lat-th1;
            [DTH R] = meshgrid(dth,axial);
            t2 = t([1:5 8:end]);
            T = repmat(reshape(t2(:),1,1,[]),size(R,1),size(R,2));
            DX = R.*sind(DTH) - repmat(A.*cosd(th1),size(R,1),1).*sind(DTH);
            DX3 = repmat(DX,[1 1 size(arfidata2,3)-1]);
            V3 = diff(arfidata2,1,3)./diff(T,1,3);
            T3 = 0.5*T(:,:,1:end-1)+0.5*T(:,:,2:end);
            t3 = reshape(T3(1,1,:),[],1);
            DXDT = abs(DX3)./T3;
            fs = 1/(mean(diff(t*1e-3)));
            V4 = reshape(V3,size(V3,1),par.nBeams,[]);
            V4 = V4 - repmat(mean(V4,2),[1 size(V4,2) 1]);
            V4 = reshape(V4,size(V4,1),par.nBeams*par.numBeamGroups*par.numAcq,[]);
            passbandHz =[25 1000];
            [BB,AA] = butter(3, min(max(passbandHz/(fs/2),1e-12),1));
            V5 = flipdim(filter(BB,AA,flipdim(filter(BB,AA,V4,[],3),3),[],3),3);
            VTEST = convn(V5,ones(5,1)/5,'same');
            %VTEST((T3>1)&(abs(DX3)./(T3-1))<sws.dxdtrange(1)) = nan;
            VTEST(DXDT<sws.dxdtrange(1)) = nan;
            VTEST(DXDT>sws.dxdtrange(2)) = nan;
            xrange = [0.5 3];
            VTEST(abs(DX3)<=xrange(1)) = nan;
            VTEST(abs(DX3)>=xrange(2)) = nan;
            VTEST(T<=0) = nan;
            [pk idx] = nanmax(VTEST,[],3);
            idx2 = mediannan(idx,[3 3]);
            DT = t3(round(idx2));
            DT(isnan(pk)) = nan;
            
            DT = reshape(DT,size(VTEST,1),par.nBeams,[]);
            for i = 1:size(DT,3);
                  DT(:,:,i) = mediannan(DT(:,:,i),[3 3]);
            end
            DX = reshape(DX,size(VTEST,1),par.nBeams,[]);
            M = ones(size(DX));
            M(abs(DX)<=xrange(1)) = nan;
            M(abs(DX)>=xrange(2)) = nan;
            M(isnan(DT)) = nan;
            
            DTDX = (nanmean(DT.*abs(DX).*M,2)-(nanmean(DT.*M,2).*nanmean(abs(DX).*M,2)))./(nanmean(M.*abs(DX).^2,2)-nanmean(abs(DX).*M,2).^2);
            T0 = nanmean(DT.*M,2)-DTDX.*nanmean(abs(DX).*M,2);
            
            CT = abs(DX)./(DT-repmat(min(2,max(0,T0)),[1 par.nBeams 1]));
            
            %CT = abs(DX)./DT;
            CT = reshape(CT,size(CT,1),par.nBeams,[]);
            CT(CT<sws.dxdtrange(1)) = nan;
            CT(CT>sws.dxdtrange(2)) = nan;
            figure(figHandles.SWEI);clf;
            %im2 = surf(X(:,1:par.nBeams),Y(:,1:par.nBeams),0*R(:,1:32),nanmean(CT,3));set(im2,'edgecolor','none');axis ij;axis image;view([0 90]);
            im2 = imsurf(lat(1,1:par.nBeams),axial,A,nanmedian(CT,3));
            axis(rt.axisLimits);
            caxis([rt.dispLim])
            cb = colorbar('location','northoutside');
            xlabel(cb,'SWS (m/s)')
            xlabel('x (mm)');
            ylabel('z (mm)');
            set(gca,'color','k')
            set(figHandles.SWEI,'Name','SWS','NumberTitle','off')
            grid off
            box on
            hold on
            bim3 = image(blat,baxial,bimg(:,:,[1 1 1]));
            set(bim3,'alphadata',min(0.8,1-bimg(:,:,1).^0.5));
            set(bim3,'alphadata',min(0.9,bimg(:,:,1)<20/255))
            if ~rt.bMask;set(bim3,'visible','off');end
            title('')
            caxis([0 sws.dxdtrange(2)])
            
            figure(110)
            [ii jj] = sort(par.pushbeamNum);
            im3 = imsurf(lat(1,1:par.nBeams),axial,A,squeeze(DTDX(:,jj)),[0 1]);
            set(110,'Name','DTDX','NumberTitle','off');
            cb = colorbar('location','northoutside');
            set(gca,'color','k')
            
            diffn = 1;
            CTSMURF = squeeze(nanmedian(diff(DX(:,:,jj),diffn,3)./diff(convn(DT(:,:,jj),ones(1,1,3)/3,'same').*sign(DX(:,:,jj)),diffn,3),2));
            figure(111)
            [ii jj] = sort(par.pushbeamNum);
            thsmurf = filter2(ones(1,diffn+1)/(diffn+1),lat(1,1:par.nBeams),'valid');
            im4 = imsurf(thsmurf,axial,A,CTSMURF,[0 sws.dxdtrange(2)]);
            set(111,'Name','SMURF','NumberTitle','off');
            cb = colorbar('location','northoutside');
            set(gca,'color','k')

            
        end
        
end
