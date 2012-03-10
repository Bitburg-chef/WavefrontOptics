function wvf = wvfSet(wvf,parm,val,varargin)
%Get wavefront structure parameters and derived properties
%
%     wvf = wvfSet(wvf,parm,val,varargin)
%
% Wavefront properties are either stored as parameters or computed from those
% parameters. We generally store only unique values and  calculate all
% derived values.
%
% A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used, where
% unit specifies the spatial scale of the returned value:  'm', 'cm', 'mm',
% 'um', 'nm'.  Default is always meters ('m').  (REALLY?  NOT YET).
%
% Examples:
%
%
% Parameters
%  % General
%     'name'
%     'type'
%     'pupilsize'
% 
%  % Spectral matters
%     'wave'
%     'infocuswavelength'
% 
%  % Pupil parameters
%     'calculatedpupil'
%     'measuredpupil'
%     
%  % Field?
%     'fieldsizepixels'
%     'fieldsizemm'
%
%  % Focus parameters
%     'zcoef'
%     'defocusdiopters'
%     'weightspectrum'
%     'sceparams'
%     'psf'
%
% (c) Wavefront Toolbox Team 2011

if ~exist('parm','var') || isempty(parm), error('Parameter must be defined.'); end
if ~exist('val','var'), error('val must be defined.'); end

parm = ieParamFormat(parm);

switch parm
    case 'name'
        wvf.name = val;
    case 'type'
        wvf.type = val;
    case 'pupilsize'

        % Spectral matters
    case {'wave','wavelist','wavelength'}
        % Vector of wavelengths  e.g., w = StoWls([400 5 61])
        wvf.wls = val(:);
    case {'infocuswavelength','nominalfocuswl'}
        % In focus wavelength nm.  Single value.
        wvf.nominalFocusWl = val; 

        % Pupil parameters
    case {'calculatedpupil','calculatedpupildiameter'}
        % Pupil diameter in mm - must be smaller than measurements
        wvf.calcpupilMM = val;               
    case {'measuredpupil','measuredpupildiameter'}
        % Largest measured pupil diameter in mm
        wvf.measpupilMM = val;               
    case {'pupilfunction','pupilfunc'}
        wvf.pupilfunc = val;
        
        % Field?
    case 'fieldsizepixels'
        wvf.sizeOfFieldPixels = val;          % In pixels?  No units?
    case 'fieldsizemm'
        wvf.sizeOfFieldMM = val;
        
        % Focus parameters
    case {'zcoef','zcoeffs'}
        wvf.zcoeffs = val;
    case 'defocusdiopters'
        % Stored this way, but can be accessed in distance units.
        wvf.defocusDiopters = val;           % Defocus
    case 'strehl'
        % Not sure why this is set.  It is derived
        wvf.strehl = val;
        
    case 'defocusmicrons'
        % Hmmm.  Someone decided to have two ways of specifying defocus.
        % This is unfortunate.  Anyway, we are supposed to be able to set
        % the defocus in microns as well as in diopters.  This happens in
        % wvfComputePSF.
        wvf.defocusMicrons = val;
    case 'weightspectrum'
        wvf.weightingSpectrum = val;         % Defocus

        % Special cases
    case {'sceparams','stilescrawford'}
        % The structure of sce should be:
        %
        wvf.sceParams = val;
        
        % Derived parameters
    case 'psf'
        wvf.psf = val;

        % These pixel related measures computed in wvfComputePSF.  Not sure
        % what they are.
    case {'areapix'}
        % Don't know what this is.
        wvf.areapix = val;
    case {'areapixapod'}
        % Don't know what this is.
        wvf.areapixapod = val; 
        
    otherwise
        error('Unknown parameter %s\n',parm);
end

return
