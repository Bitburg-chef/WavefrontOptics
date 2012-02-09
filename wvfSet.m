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
    case {'wave','wavelist'}
        wvf.wls = val;
    case 'infocuswavelength'
        wvf.nominalFocusWl = val;            % In focus wavelength (nm)

        % Pupil parameters
    case 'calculatedpupil'
        wvf.calcpupilMM = val;               % Default pupil??? radius?
    case {'measuredpupil','measuredpupildiameter'}
        wvf.measpupilMM = val;               % Default pupil diameter?
    
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

    case 'weightspectrum'
        wvf.weightingSpectrum = val;         % Defocus

        % Special cases
    case 'sceparams'
        % The structure of sce should be:
        %
        wvf.sce = val;
        
        % Derived parameters
    case 'psf'
        wvf.psf = val;

    otherwise
        error('Unknown parameter %s\n',parm);
end

return
