function [I,Q] = genHarmonicSummedData(I,Q,par)

posPulse = nan(1,par.ensemble);
negPulse = nan(1,par.ensemble);

% account for references
for i = 1:par.nref-1
    posPulse(i) = i+mod(i-1,2);
    negPulse(i) = i+mod(i,2);
end
posPulse(par.nref) = posPulse(par.nref-1);
negPulse(par.nref) = negPulse(par.nref-1);

% account for pushes within the ensemble (sum same time step)
if par.separateFocalZoneAcqs
    for i = (par.nref+1):(par.npush+par.nref)
        posPulse(i) = i;
        negPulse(i) = i;
    end
else
    for i = (par.nref+1):(par.npush*length(par.pushFocalDepth)+par.nref)
        posPulse(i) = i;
        negPulse(i) = i;
    end
end

% account for tracks
for i = (i+1):par.ensemble
    posPulse(i) = i+mod(i-1,2);
    negPulse(i) = i+mod(i,2);
end

% account for anything that was greater than the ensemble length
posPulse(posPulse>par.ensemble) = par.ensemble-mod(par.ensemble-1,2);
negPulse(negPulse>par.ensemble) = par.ensemble-mod(par.ensemble,2);

I = (I(:,:,posPulse)+I(:,:,negPulse))./2;
Q = (Q(:,:,posPulse)+Q(:,:,negPulse))./2;

end