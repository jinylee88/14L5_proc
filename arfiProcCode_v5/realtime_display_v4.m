
figure(figHandles.IQBMode);clf
set(figHandles.IQBMode,'Name','IQ B-Mode','NumberTitle','off')
ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',figHandles.IQBMode);
%image(blat,baxial,bimg(:,:,[1 1 1]));
imsurf(Th1,axial,apex,bmodedata0(:,:,[1 1 1]));
xlabel('x (mm)')
ylabel('z (mm)')
%title(timeStamp)
colormap gray
colorbar('location','northoutside')
axis image
axis(rt.axisLimits);

switch swifParams.imagingModeShort
    case 'arfi'
        figure(figHandles.ARFI);clf
        ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',figHandles.ARFI);
        im0 = imagesc([-0.1 0.1],[0 0.1],arfidata(:,:,12));
        set(im0,'visible','off');
        hold on
        tidx = find(t<=rt.t_arfi_ms,1,'last');
        if rt.depthCorrect
        R0 = sqrt((Z-6).^2+X.^2);
        normalization = min(2,max(1,(R0*0.5).^0.5));
        im  = imsurf(theta,axial,apex,double(arfidata(:,:,tidx).*normalization));
        else
        im = imsurf(theta,axial,apex,double(arfidata(:,:,tidx)));
        end
        set(ax,'color','k')
        xlabel('x (mm)')
        ylabel('z (mm)')
        cb = colorbar('NorthOutside');
        set(cb,'XColor','k')
        if rt.depthCorrect
        xlabel(cb,'Normalized Displacement (\mum)','color','k')
        else
        xlabel(cb,'Displacement (\mum)','color','k')
        end
        bim = image(blat,baxial,bimg(:,:,[1 1 1]));
        bim2 = image(blat,baxial,bimg(:,:,[1 1 1]));
        set(ax,'Children',[bim2 im bim im0])
        if ~rt.bMask;
            set(bim2,'visible','off');
        else
            set(bim2,'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8)))
        end
        axis image
        axis(rt.axisLimits)
        caxis(rt.dispLim)
        set(figHandles.ARFI,'Name','ARFI','NumberTitle','off');
        drawnow
        
        if flags.calcSWS && strcmpi(par.mode,'mmode');
            for i = 1:length(swsmethod)
            figure(figHandles.SWEI(i));clf;
            ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',figHandles.SWEI(i));
            eval(sprintf('CT = CT_%s;',upper(swsmethod{i})));
            im2 = imsurf(theta,axial,apex,nanmedian(CT,3));
            axis(rt.axisLimits);
            set(gca,'color','k')
            caxis([rt.dispLim])
            cb = colorbar('location','northoutside');
            set(cb,'xColor','k')
            xlabel(cb,'SWS (m/s)','color','k')
            xlabel('x (mm)');
            ylabel('z (mm)');
            set(figHandles.SWEI(i),'Name',['SWS_' upper(swsmethod{i})],'NumberTitle','off')
            grid off
            box on
            hold on
            if rt.bMask
            bim3 = image(blat,baxial,bimg(:,:,[1 1 1]));
            set(bim3,'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8)));
            end
            title('')
            caxis([0 sws.dxdtrange(2)])
            end
        end
    case 'mmode'
            
            if flags.displayflag>1
            for i = 1:length(swsmethod)
            figure(figHandles.SWEI(i));clf;
            ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',figHandles.SWEI(i));
            eval(sprintf('CT = CT_%s;',upper(swsmethod{i})));
            im2 = imsurf(theta,axial,apex,nanmedian(CT,3));
            axis(rt.axisLimits);
            set(gca,'color','k')
            caxis([rt.dispLim])
            cb = colorbar('location','northoutside');
            set(cb,'xColor','k')
            xlabel(cb,'SWS (m/s)','color','k')
            xlabel('x (mm)');
            ylabel('z (mm)');
            set(figHandles.SWEI(i),'Name',['SWS_' upper(swsmethod{i})],'NumberTitle','off')
            grid off
            box on
            hold on
            bim3 = image(blat,baxial,bimg(:,:,[1 1 1]));
            set(bim3,'alphadata',min(0.8,1-bimg(:,:,1).^0.5));
            set(bim3,'alphadata',min(0.9,bimg(:,:,1)<20/255))
            if ~rt.bMask;set(bim3,'visible','off');end
            title('')
            caxis([0 sws.dxdtrange(2)])
            end
            end
           
            if frmii == 1
                figure(figHandles.MMODE);
                clf(figHandles.MMODE);
                ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',figHandles.MMODE);
                im4 = imagesc(Frameindices,axial,zeros(length(axial),length(Frameindices)));
                set(ax,'ylim',rt.axisLimits(3:4));
                set(ax,'color','k')
                cb = colorbar('location','northoutside');
                set(cb,'xColor','k')
                xlabel(cb,'SWS (m/s)','color','k')
                xlabel('Frame #');
                ylabel('z (mm)');
                set(figHandles.MMODE,'Name','SWS-MMODE','NumberTitle','off')
                grid off
                box on
                caxis([0 sws.dxdtrange(2)])
                hold on
                bim4 = image(Frameindices,baxial,zeros(length(axial),length(Frameindices)));
                if ~rt.bMask;set(bim4,'visible','off');end
            else
                mmode = zeros(length(axial),length(Frameindices));
                mmode(:,1:frmii) = squeeze(mean(Bmodedata0(:,round(size(Bmodedata0,2)/2)+[-2:2],:),2));
                cdata = zeros(length(axial),length(Frameindices));
                cdata(:,1:frmii) = squeeze(nanmedian(SWS,2));
                set(bim4,'cdata',mmode,'alphadata',min(0.9,mmode<20/255))
                set(im4,'cdata',cdata,'alphadata',~isnan(cdata));
                drawnow;
            end
            

            
            
end
