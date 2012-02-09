function val = wvfGet(wvf,parm,varargin)
%Get wavefront structure parameters and derived properties
%
%     val = wvfGet(wvf,parm,varargin)
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
    case 'nwave'
        val = length(wvf.wls);
    case 'infocuswavelength'
        wvf.nominalFocusWl = 650;            % In focus wavelength (nm)

        % Pupil parameters
    case 'calculatedpupil'
        val = wvf.calcpupilMM;               % Default pupil??? radius?
    case 'measuredpupil'
        % Default is in millimeters
        % wvfGet(wvf,'measured pupil','mm');
        val = wvf.measpupilMM;               % Default pupil diameter?
        % Scale for unit
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val/1000)*ieUnitScaleFactor(varargin{1});
        end
   
        % Field?
    case {'fieldsizepixels','npixels','fieldsize'}
        % In pixels?  No units?  Why not a distance or an angle or
        % something?
        val = wvf.sizeOfFieldPixels;         
    case 'fieldsizemm'
        val = wvf.sizeOfFieldMM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val/1000)*ieUnitScaleFactor(varargin{1});
        end
        
        % Focus parameters
    case {'zcoeffs'}
        % wvfGet(wvf,'zcoef',list)
        % wvfGet(wvf,'zcoef',4)
        if isempty(varargin)
            val = wvf.zcoeffs;
        else
            val = wvf.zcoeffs(varargin{1});
        end
        
    case 'defocusdiopters'
        val = wvf.defocusDiopters;           % Defocus
    case {'defocusmicrons','defocusdistance'}
        % The defocus in distance rather than diopters
        % The default is microns.
        % wvfGet(wvfP,'defocus distance','mm');
        val = wvfGetDefocusFromWavelengthDifference(wvf);
        if isempty(varargin), return
        else 
            % Convert to meters and then scale
            val = (val/10^6)*ieUnitScaleFactor(varargin{1});
        end
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
