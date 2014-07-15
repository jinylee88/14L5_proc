function par = checkParams(par)

% check existence of params
fieldName = {'kernelLength','interpFactor','pushFocalDepth','pushFreq','pushFnum','txTypeIndex','lambda','pushPRF','ref_idx'};
for i = 1:length(fieldName)
    if ~isfield(par,fieldName{i}),error('Could not find %s in par\n',fieldName{i});end
end

% add unit information
par.units.system = 'MKS';
par.units.t = 'ms';
par.units.pushFocalDepth = 'mm';
par.units.trackFocalDepth = 'mm';
par.units.lambda = 'mm';
par.units.kernelLength = 'wavelength';
par.units.pushPRF = 'Hz';

if strcmp(par.probeType,'linear')
    par.units.axial = 'mm';
    par.units.lat = 'mm';
elseif max(strcmp(par.probeType,{'curvilinear','phased'}))
    par.units.radial = 'mm';
    par.units.angular = 'rad';
else
    error('Unknown probe type')
end



end