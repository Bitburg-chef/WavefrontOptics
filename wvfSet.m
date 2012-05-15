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
% History:
%   4/29/12  dhb  Allow wls as a synonym for wavelength, because that was
%                 my first guess given the name of the field.
%
% (c) Wavefront Toolbox Team 2011, 2012

% Programming Notes
%   The Strehl ratio, http://en.wikipedia.org/wiki/Strehl_ratio
%

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
    case {'wave','wavelist','wavelength','wls'}
        % Normally just a vector of wavelengths
        % but allow SToWls case for DHB.
        % wvfSet(wvfP,'wave',400:10:700)  OR
        % wvfSet(wvfP,'wave',[400 10 31])
        if size(val,2) == 3 && size(val,1) == 1 % SToWls case
            % Row vector with 3 entries.
            wls = SToWls(val);
            wvf.wls = wls;
        else  % A column vector case
            wvf.wls = val(:);
        end
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
        % wvfSet(wvf,'pupil function',pf) - pf is a cell array of pupil
        % functions, one for each wavelength
        %
        % pupil function is a point spread function matrix for the
        % wave(idx). wvfSet(wvf,'pupil function',pf,[idx])
        %
        % Convert between pupil function and psf using wvfComputePSF.

        if isempty(varargin)  % Cell array of pupilfuncs
            % This is a cell array of pupil functions if there are multiple
            % wavelengths, or just a matrix
            if iscell(val)
                % Check cell array dimension
                n = length(val); nWave = wvfGet(wvf,'nWave');
                if n ~= nWave
                    error('pupilfunc dim (%d) ~= nWave (%d)', n, nWave);            
                end
                wvf.pupilfunc = val;
            else  % Just a matrix, not a cell array.  
                % No idx specified, so we put it in the first cell.
                warning('WVFSET:pupilfuncset','Assigning pupil function to first cell array dim');
                wvf.pupilfunc{1} = val;
            end
        else  % The wavelength index was sent in
            % wvfSet(wvf,'pupilfunc',val,idx)
            % This is the pupilfunc for wave(idx)
            idx = varargin{1};
            nWave = wvfGet(wvf,'n wave');
            if idx > nWave, error('idx (%d) > nWave (%d)',idx,nWave);
            else            wvf.pupilfunc{idx} = val;
            end
        end
        
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
        % Zernicke coefficients.  This should be 65 terms, and we
        % over-write the number that are in val.
        wvf.zcoeffs(1:length(val)) = val;
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
        
        % PSF parameters
    case 'psf'
        % wvfSet(wvf,'psf',psf) - psf is a cell array of point spreads, one
        % for each wavelength
        %
        % psf is a point spread function matrix for the wave(idx).
        % wvfSet(wvf,'psf',psf,[idx])
        %
        % Point spread function.  Handled as cell array because of
        % potential differences in dimension due to wavelength.  See
        % pupilfunc.  
        if isempty(varargin)  % Cell array of pupilfuncs
            % This is a cell array of pupil functions if there are multiple
            % wavelengths, or just a matrix
            if iscell(val)
                % Check cell array dimension
                n = length(val); nWave = wvfGet(wvf,'nWave');
                if n ~= nWave, error('psf dim (%d) ~= nWave (%d)', n, nWave);
                else           wvf.psf = val;
                end
            else  % Just a matrix, not a cell array.  
                % No idx specified, so we put it in the first cell.
                warning('WVFSET:pupilfuncset','Assigning pupil function to first cell array dim');
                wvf.psf{1} = val;
            end
        else  % The wavelength index was sent in
            % wvfSet(wvf,'pupilfunc',val,idx)
            % This is the pupilfunc for wave(idx)
            idx = varargin{1};
            nWave = wvfGet(wvf,'n wave');
            if idx > nWave, error('idx (%d) > nWave (%d)',idx,nWave);
            else            wvf.psf{idx} = val;
            end
        end
        
%         % These pixel related measures computed in wvfComputePupilFunction
%         % and stored in wvComputePSF.  Not sure what they are and whether
%         % they need to be here.  Perhaps we could move them from that
%         % function into here?  They probably should never be set, but they
%         % should always be a get based on pupilfunc.
%         % BW, May 2012.
%     case {'areapix'}
%         % Not sure about the physical significance of this If we know the
%         % area of each pixel, we can use this to calculate the area covered
%         % by the pupil function. It is computed for the first time in
%         % wvfComputePupilFunction as numel(pupilfunc))
%         wvf.areapix = val;
%     case {'areapixapod'}
%         % Not sure about the physical significance of this Something like
%         % the area underneath the absolute value of the pupil function It
%         % is a vector the same length as wavelength. It is computed for the
%         % first time in wvfComputePupilFunction sum(sum(abs(pupilfunc)))
%         wvf.areapixapod = val;
    otherwise
        error('Unknown parameter %s\n',parm);
end

return
