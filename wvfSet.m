function wvf = wvfSet(wvf,parm,val,varargin)
%Get wavefront structure parameters and derived properties
%
%     wvf = wvfSet(wvf,parm,val,varargin)
%
% Wavefront properties are either stored as parameters or computed from
% those parameters. We generally store only unique values and  calculate
% all derived values.
%
% A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used; the
% parameter 'unit' specifies the spatial scale of the returned value:  'm',
% 'cm', 'mm', 'um', 'nm'.  There is not always a default in meters ('m'),
% as there should be.
%
% Parameter names can be written with spaces and upper/lower case.  The
% strings are converted to lower case and all the spaces are removed by
% this routine.
%
% Examples:
%   wvf = wvfSet(wvf,'measured pupil',3);   % Default is mm, I think
%   wvf = wvfSet(wvf,'stiles crawford',sce);
%   wvf = wvfSet(wvf,'psf',psf);
%   wvf = wvfSet(wvf,'infocus wavelength',550);
%
% Parameters
%  % General
%     'name' - Name of this object
%     'type' - 
% 
%  % Spectral matters
%     'wave'
%     'infocuswavelength'
% 
%  % Pupil parameters
%     'calculatedpupil'
%     'measuredpupil'
%     'fieldSampleSize'
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
        % This specific object name
        wvf.name = val;
    case 'type'
        % Should always be wvf.  
        wvf.type = val;

        % Spectral matters
    case {'wave','wavelist','wavelength'}
        % Vector of wavelengths  e.g., w = SToWls([400 5 61])
        % or single wavelength e.g. w = SToWls([550 1 1]) or w = 550;
        wls = SToWls(val);
        wvf.wls = wls;
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
        % pupil field amplitude/phase function
        wvf.pupilfunc = val;
    case {'fieldsamplesize','fieldsamplesizemmperpixel'}
        % pixel size of sample pupil field, used to calculate number of
        % pixels of computed field
        wvf.fieldSampleSizeMMperPixel = val;         
        % need to make sure field size is integer multiple of sample size
        % (so that there are an integer number of pixels)
        nPixels = ceil(wvf.sizeOfFieldMM/val);
        wvf.sizeOfFieldMM = val*nPixels;
    case 'fieldsizemm'
        % total size of computed field in pupil plane, used to calculate
        % number of pixels of computed field
        wvf.sizeOfFieldMM = val;
        % need to make sure field size is integer multiple of sample size
        % (so that there are an integer number of pixels)
        nPixels = ceil(val/wvf.fieldSampleSizeMMperPixel);
        wvf.fieldSampleSizeMMperPixel = val/nPixels;
        
        % Focus parameters
    case {'zcoef','zcoeffs'}
        % Zernicke coefficients.
        wvf.zcoeffs = val;
    case 'defocusdiopters'
        % Does not look like defocus is ever stored in diopters in other
        % functions. Only for user to set defocus diopters as an
        % alternative to using 4th term of zernike coeffs.
        wvf.defocusDiopters = val;           % Defocus
    case 'strehl'
        % Not sure why this is set.  It is derived
        wvf.strehl = val;
        
    case 'defocusmicrons'
        % Hmmm.  Someone decided to have two ways of specifying defocus.
        % This is unfortunate.  Anyway, we are supposed to be able to set
        % the defocus in microns as well as in diopters.  This happens in
        % wvfComputePSF.
        % defocusMicrons is set by wvfGetDefocusFromWavelengthDifference.
        % It calculates additional longitudinal chromatic aberration (LCA)
        % from using non-nominal focus wavelengths and adds it to 
        % user-specified defocus diopters, then converts it all to um.
        % wvfComputePSF adds in this defocus microns (same units as rest of
        % pupil function calculations) to zcoeff(4)/defocus if present.
        % KP 3/12/12
        wvf.defocusMicrons = val;
    case 'weightspectrum'
        % Used for calculating defocus to white light with many spectral
        % terms. 
        wvf.weightingSpectrum = val;         % Defocus

        % Special cases
    case {'sceparams','stilescrawford'}
        % The structure of sce is defined in sceCreate
        wvf.sceParams = val;
        
        % Derived parameters
    case 'psf'
        % Point spread function.  Not sure how many different wavelengths
        % are handled.  Define, please.
        wvf.psf = val;

        % These pixel related measures computed in wvfComputePSF.  Not sure
        % what they are and whether they need to be here.  Perhaps we could
        % move them from that function into here?
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
