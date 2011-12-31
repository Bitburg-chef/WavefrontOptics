function val = wvfGet(wvf,parm,varargin)
%Get wavefront structure parameters and derived properties
%
%     val = wvfGet(scene,parm,varargin)
%
% Wavefront properties are either stored as parameters or computed from those
% parameters. We generally store only unique values and  calculate all
%  derived values.
%
%  A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used, where
%  unit specifies the spatial scale of the returned value:  'm', 'cm', 'mm',
%  'um', 'nm'.  Default is always meters ('m').  (REALLY?  NOT YET).
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

% Default is empty when the parameter is not yet defined.
val = [];

parm = ieParamFormat(parm);

switch parm
    case 'name'
        val = wvf.name;
    case 'type'
        val = wvf.type;
    case 'pupilsize'

        % Spectral matters
    case 'wave'
        val = wvf.wls;
    case 'infocuswavelength'
        wvf.nominalFocusWl = 650;            % In focus wavelength (nm)

        % Pupil parameters
    case 'calculatedpupil'
        val = wvf.calcpupilMM;               % Default pupil??? radius?
    case 'measuredpupil'
        val = wvf.measpupilMM;               % Default pupil diameter?
    
        % Field?
    case 'fieldsizepixels'
        val = wvfP.sizeOfFieldPixels;          % In pixels?  No units?
    case 'fieldsizemm'
        val = wvfP.sizeOfFieldMM;
        
        % Focus parameters
    case 'zcoef'
        val = wvf.zcoeffs;
    case 'defocusdiopters'
        val = wvf.defocusDiopters;           % Defocus
    case 'weightspectrum'
        val = wvf.weightingSpectrum;         % Defocus

        % Special cases
    case 'sceparams'
        val = wvf.sce;
        
        % Derived parameters
    case 'psf'
        val = wvf.psf;

    otherwise
        error('Unknown parameter %s\n',parm);
end

return
