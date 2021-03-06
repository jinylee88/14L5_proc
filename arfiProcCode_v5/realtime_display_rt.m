if fileidx == 1
figure(options.display.figureHandles.IQBMode);clf
set(options.display.figureHandles.IQBMode,'Name','IQ B-Mode','NumberTitle','off')
ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.IQBMode);
%image(blat,baxial,bimg(:,:,[1 1 1]));
switch par.probe
    case 'linear'
    BMODE_handle = image(blat,axial,bmodedata0(:,:,[1 1 1]));
    case 'phased'
    BMODE_handle = imsurf(blat,axial,apex,bmodedata0(:,:,[1 1 1]));
end
xlabel('x (mm)')
ylabel('z (mm)')
%title(timeStamp)
colormap gray
colorbar('location','northoutside')
axis image
axis(options.display.axisLimits);
else
    set(BMODE_handle,'cdata',bmodedata0(:,:,[1 1 1]));
    drawnow
end


switch swifParams.imagingModeShort
    case 'arfi'
        figure(options.display.figureHandles.ARFI);clf
        ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.ARFI);
        im0 = imagesc([-0.1 0.1],[0 0.1],arfidata(:,:,12));
        set(im0,'visible','off');
        hold on
        tidx = find(t<=options.display.t_arfi_ms,1,'last');
        if options.display.depthCorrect
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
        if options.display.depthCorrect
        xlabel(cb,'Normalized Displacement (\mum)','color','k')
        else
        xlabel(cb,'Displacement (\mum)','color','k')
        end
        bim = image(blat,baxial,bimg(:,:,[1 1 1]));
        bim2 = image(blat,baxial,bimg(:,:,[1 1 1]));
        set(ax,'Children',[bim2 im bim im0])
        if ~options.display.bmodeMask;
            set(bim2,'visible','off');
        else
            set(bim2,'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8)))
        end
        axis image
        axis(options.display.axisLimits)
        caxis(options.display.dispLim)
        set(options.display.figureHandles.ARFI,'Name','ARFI','NumberTitle','off');
        drawnow
        
        if flags.calcSWS && strcmpi(par.mode,'mmode');
            for i = 1:length(swsmethod)
            figure(options.display.figureHandles.SWEI(i));clf;
            ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.SWEI(i));
            eval(sprintf('CT = CT_%s;',upper(swsmethod{i})));
            im2 = imsurf(theta,axial,apex,nanmedian(CT,3));
            axis(options.display.axisLimits);
            set(gca,'color','k')
            caxis([options.display.arfiLimit])
            cb = colorbar('location','northoutside');
            set(cb,'xColor','k')
            xlabel(cb,'SWS (m/s)','color','k')
            xlabel('x (mm)');
            ylabel('z (mm)');
            set(options.display.figureHandles.SWEI(i),'Name',['SWS_' upper(swsmethod{i})],'NumberTitle','off')
            grid off
            box on
            hold on
            if options.display.bmodeMask
            bim3 = image(blat,baxial,bimg(:,:,[1 1 1]));
            set(bim3,'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8)));
            end
            title('')
            caxis([options.display.dxdtLimit])
            end
        end
    case 'mmode'
            
            if flags.displayflag>1
            for i = 1:length(swsmethod)
            figure(options.display.figureHandles.SWEI(i));clf;
            ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.SWEI(i));
            eval(sprintf('CT = CT_%s;',upper(swsmethod{i})));
            im2 = imsurf(theta,axial,apex,nanmedian(CT,3));
            axis(options.display.axisLimits);
            set(gca,'color','k')
            caxis([options.display.arfiLimit])
            cb = colorbar('location','northoutside');
            set(cb,'xColor','k')
            xlabel(cb,'SWS (m/s)','color','k')
            xlabel('x (mm)');
            ylabel('z (mm)');
            set(options.display.figureHandles.SWEI(i),'Name',['SWS_' upper(swsmethod{i})],'NumberTitle','off')
            grid off
            box on
            hold on
            bim3 = image(blat,baxial,bimg(:,:,[1 1 1]));
            set(bim3,'alphadata',min(0.8,1-bimg(:,:,1).^0.5));
            set(bim3,'alphadata',min(0.9,bimg(:,:,1)<20/255))
            if ~options.display.bmodeMask;set(bim3,'visible','off');end
            title('')
            caxis(options.display.dxdtLimit)
            end
            end
           
            if frmii == 1
                figure(options.display.figureHandles.MMODE);
                clf(options.display.figureHandles.MMODE);
                ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.MMODE);
                im4 = imagesc(Frameindices,axial,zeros(length(axial),length(Frameindices)));
                set(ax,'ylim',options.display.axisLimits(3:4));
                set(ax,'color','k')
                cb = colorbar('location','northoutside');
                set(cb,'xColor','k')
                xlabel(cb,'SWS (m/s)','color','k')
                xlabel('Frame #');
                ylabel('z (mm)');
                set(options.display.figureHandles.MMODE,'Name','SWS-MMODE','NumberTitle','off')
                grid off
                box on
                caxis(options.display.dxdtLimit)
                hold on
                bim4 = image(Frameindices,baxial,zeros(length(axial),length(Frameindices)));
                if ~options.display.bmodeMask;set(bim4,'visible','off');end
            else
                mmode = zeros(length(axial),length(Frameindices));
                mmode(:,1:frmii) = squeeze(mean(Bmodedata0(:,round(size(Bmodedata0,2)/2)+[-2:2],:),2));
                cdata = zeros(length(axial),length(Frameindices));
                cdata(:,1:frmii) = squeeze(nanmedian(SWS,2));
                set(bim4,'cdata',mmode,'alphadata',min(0.9,mmode<20/255))
                set(im4,'cdata',cdata,'alphadata',~isnan(cdata));
                drawnow;
            end
    case 'arfiswei'
        
        switch par.probe
            case 'linear'
                if fileidx == 1
              figure(options.display.figureHandles.ARFI);clf
        ax = subplot('position',[0.075 0.05 0.85 0.9],'parent',options.display.figureHandles.ARFI);
        tidx = find(t<=options.display.arfiFrameTms,1,'last');
        im = imagesc(mean(lat,2),axial,arfidata(:,:,tidx));
        hold on
        set(ax,'color','k')
        xlabel('x (mm)')
        ylabel('z (mm)')
        cb = colorbar('NorthOutside');
        set(cb,'XColor','k')
        xlabel(cb,'Displacement (\mum)','color','k')
        bim = image(blat,axial,bmodedata0(:,:,[1 1 1]));
        bim2 = image(blat,axial,bmodedata0(:,:,[1 1 1]));
        set(ax,'Children',[bim2 im bim])
        if ~options.display.bmodeMask;
            set(bim2,'visible','off');
        else
            set(bim2,'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8)))
        end
        axis image
        axis(options.display.axisLimits)
        caxis(options.display.arfiLimit)
        set(options.display.figureHandles.ARFI,'Name','ARFI','NumberTitle','off');
        drawnow
                else
                    set(im,'cdata',arfidata(:,:,tidx))
                    set(bim,'cdata',bmodedata0(:,:,[1 1 1]));
                   if options.display.bmodeMask 
                       set(bim2,'cdata',bmodedata0(:,:,[1 1 1],'alphadata',imclose(imopen(min(0.9,bimg(:,:,1)<20/255),strel('disk',8)),strel('disk',8))));
                   end
                   drawnow
                end
        end
            
end
