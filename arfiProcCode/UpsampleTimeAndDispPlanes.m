
%==========================================================================

function [planesup,tup]=UpsampleTimeAndDispPlanes(planes,t,desiredPRF)

        %
        % upsample t and dispPlanes to give tup and dispPlanesup.
        %   Input planes can be 2D or 3D.  desiredPRF in kHz, t in ms.
        %

    expPRF = 1/mean(diff(t));
    upsamp = round(desiredPRF/expPRF);
    
    deltat = mean(diff(t))/upsamp;
    tup=t(1):deltat:max(t);

    ndims=length(size(planes));
    
    switch ndims
        case 2
            planesup=UpsampleDispPlanes2D(planes,t,tup);
        case 3
            planesup=UpsampleDispPlanes3D(planes,t,tup);
        otherwise
            error('planes must be 2D or 3D')
    end
end
    
%==========================================================================    

function planesup=UpsampleDispPlanes2D(planes,t,tup)

    [nlats,ntimes]=size(planes);
    planesup=zeros(nlats,length(tup));
    
    for ilat=1:nlats
        planesup(ilat,:)=interp1(t,planes(ilat,:),tup,'spline');
    end    
end

%==========================================================================

function planesup=UpsampleDispPlanes3D(planes,t,tup)

    [naxial,nlats,ntimes]=size(planes);
    planesup=zeros(naxial,nlats,length(tup));

    for iaxial=1:naxial
        for ilat=1:nlats
            temp=squeeze(planes(iaxial,ilat,:));
            planesup(iaxial,ilat,:)=interp1(t,temp,tup,'spline');
        end
    end
end

%==========================================================================
