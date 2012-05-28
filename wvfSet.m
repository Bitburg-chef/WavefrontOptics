function wvf = wvfSet(wvf,parm,val,varargin)
% wvf = wvfSet(wvf,parm,val,varargin)
%
% Get wavefront structure parameters and derived properties
%
% See also: wvfGet, wvfCreate, sceCreate, sceGet
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
%   wvf = wvfSet(wvf,'measured pupil',3);   % Default is mm
%   wvf = wvfSet(wvf,'stiles crawford',sce);
%   wvf = wvfSet(wvf,'psf',psf);
%
% Parameters:
%
%  Bookkeeping
%   'name' - Name of this object
%   'type' - Type of this object, should always be 'wvf'
%
%  Zernike coefficients and related
%   'zcoeffs' - Zernike coefficients
%   'measured pupil' - Pupil size for wavefront aberration meaurements (mm)
%   'measured wl' - Wavefront aberration measurement wavelength (nm)
%   'measured optical axis' - Measured optical axis (deg)
%   'measured observer accommodation' - Observer accommodation at aberration measurement time (diopters)
%
%  Spatial sampling parameters
%    'sample interval domain' - Which domain has sampling interval held constant with wavelength ('psf', 'pupil')
%    'spatial samples' - Number of spatial samples (pixel) for pupil function and psf
%    'ref pupil plane size' - Size of sampled pupil plane at measurement wavelength (mm)
%    'ref pupil plane sample interval' - Pixel sampling interval in pupil
%    'ref psf arcmin per sample' - Sampling interval for psf at measurment wavelength (arcminute/pixel)
%
%  Spectral
%     'calc wavelengths' - Wavelengths to compute on (nm)
%
%     'infocuswavelength'
%
%  % Pupil parameters
%     'calculatedpupil'

%
%  % Focus parameters
%     'defocusdiopters'
%     'weightspectrum'
%     'sceparams'
%     'psf'
%
% Notes:
%   5/17/12  dhb  Why is setting pupilfunc and psf allowed?
%                 This seems like it could produce all sorts of
%                 inconsistencies.  For example, various things like
%                 the strehl ratio should depend on psf, and it should
%                 be consistent with the zcoefs and pf.
%            dhb  When we pass fewer than 65 coefficients, should we zero
%                 out the higher order ones, or leave them alone?  The
%                 current code leaves them alone, which seems a little
%                 dangerous.
%            dhb  There are two underlying field sizes, one in the pupil
%                 plane and one in the plane of the retina.  The pixel
%                 dimensions are implicitly linked by the conversion between
%                 pf and psf, and our conversion code uses the same number
%                 of pixels in each representation.  One could get fancier
%                 and explicitly specify the units of each representation,
%                 and appropiately convert.  An important consideration is
%                 for the dimensions to be chosen so that both pupile
%                 function and psf are adequately sampled.
%            dhb  I tend to agree with a comment below that defocus is
%                 independently specified in too many places in the current
%                 structure.  It's possible that only the zcoef for defocus
%                 should ever change.  But before hacking it up we ought to
%                 look through the various places and ways that focus is
%                 varied and come up with a coherent first principles
%                 design.
%
% History:
%   4/29/12  dhb  Allow wls as a synonym for wavelength, because that was
%                 my first guess given the name of the field.
%   5/17/12  dhb  I agree with an earlier comment, no need to be able to set strehl.
%                 This should be a computed quantity.  I commented it out
%                 to see if anything breaks.
%            dhb  Added some checks for error conditions, modified some
%                 comments, added some notes to think about.
%
% (c) Wavefront Toolbox Team 2011, 2012

% Programming Notes
%   The Strehl ratio, http://en.wikipedia.org/wiki/Strehl_ratio

% Arg checks and parse.
%
% The switch on what we're setting is broken out into several pieces
% below to allow use of cells, and so that autoindent does something
% reasonable with our block comment style.
if ~exist('parm','var') || isempty(parm), error('Parameter must be defined'); end
if ~exist('val','var'), error('val must be defined'); end
parm = ieParamFormat(parm);

%% Initialize flags
PUPILFUNCTION_STALE = false;
PSF_STALE = false;
DIDASET = false;

%% Bookkeeping
switch parm
    % Bookkeeping
    case 'name'
        % This specific object's name
        wvf.name = val;
        DIDASET = true;
        
    case 'type'
        % Type should always be 'wvf'
        if (~strcmp(val,'wvf'))
            error('Can only set type of wvf structure to ''wvf''');
        end
        wvf.type = val;
        DIDASET = true;
end

%% Zernike coefficients and related
%
% These specify the measured (or assumed) wavefront aberrations
% in terms of a Zernike polynomial expansion.  Exanding these
% gives us the wavefront abberations in microns over the
% measured pupil.
%
% The coefficients represent measurements that were made (or
% assumed to be made) at a particular optical axis, state of
% accommodation of the observer, wavelength, and over a
% particular pupil size diameter.  We specify all of this
% size information along with the coefficients, even though
% we don't know quite how to use all of it at present.
%
% Zernike coeffs 0,1,2 (piston, verticle tilt, horizontal tilt) are
% typically 0 since they are either constant (piston) or only change the
% point spread location, not quality, as measured in wavefront aberrations.
% 65 coefficients represent the "j" single-index scheme of OSA standards
% for 10 orders of terms (each order has order+1 number of terms).
% 10 orders, including 0th order, should result in 66 total terms, but 0th
% order term, coeff(0)/piston, is neglected, resulting in 65 terms.
% Thus, sample data typically starts with coeff(3), which is +/-45 degree
% astigmatism.
switch parm
    case {'zcoeffs', 'zcoeff','zcoef'}
        % Zernicke coefficients.  This should be 65 terms, and we
        % over-write the current field with the number that are in val,
        % leaving any higher order ones alone.
        if (length(val) > 65)
            error('We do not handle more than 65 coefficients');
        end
        wvf.zcoeffs(1:length(val)) = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'measuredpupil', 'measuredpupilmm', 'measuredpupildiameter'}
        % Pupil diameter in mm over for which wavefront expansion is valid
        wvf.measpupilMM = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'measuredwl', 'measuredwavelength'}
        % Measurement wavelength (nm)
        wvf.measWlNM = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'measuredopticalaxis', 'measuredopticalaxisdeg'}
        % Measurement optical axis, degrees eccentric from fovea
        wvf.measOpticalAxisDeg = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'measuredobserveraccommodation', 'measuredobserveraccommodationdiopters'}
        % Observer accommodation, in diopters relative to relaxed state of eye
        wvf.measObserverAcommodationDiopters = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
end

%% Spatial sampling parameters
%
% In the end, we calculate using discretized sampling.  Because
% the pupil function and the psf are related by a fourier transform,
% it is natural to use the same number of spatial samples for both
% the pupil function and the corresponding psf.   The parameters
% here specify the sampling.
%
% Note that this may be done independently
% of the wavefront measurements, because the spatial scale of those
% is defined by pupil size over which the Zernike coefficients define
% the wavefront aberrations.
%
% There are a number of parameters that are yoked together here, because
% the sampling intervals in the pupil and psf domains are yoked, and because the
% overall size of the sampled quantities is determined from the sampling
% intervals and the number of pixels.
%
% As a general matter, it will drive us nuts to have the number of pixels varying with
% anything, because then we can't sum over wavelength easily.  So
% we set the number of number of pixels and size of the sampling in the
% pupil plane at the measurement wavelength, and then compute/set
% everything else as needed.
%
% Because the sampling in the pupil and psf domains is yoked, its important
% to choose values that do not produce discretization artifacts in either
% domain.  We don't have an automated way to do this, but the default
% numbers here were chosen by experience (well, really by Heidi Hofer's
% experience) to work well.  Be attentive to this aspect if you decide
% you want to change them by very much.
%
% The other thing that is tricky is that the relation between the sampling
% in the pupil and psf domains varies with wavelength. So, you can't
% easily have the sample interval stay constant over wavelength in both
% the pupil and psf domains.  You have to choose one or the other.
% We will typically force equal sampling in the psf domain, but we allow
% specification of which.
switch parm
    case {'sampleintervaldomain'}
        % Determine what's held constant with calculated wavelength.
        % Choices are 'psf' and 'pupil'
        wvf.constantSampleIntervalDomain = val;
        DIDASET = true;
        
    case {'spatialsamples', 'npixels', 'fieldsizepixels'}
        % Number of pixels that both pupil and psf planes are discretized
        % with.
        %
        % This is a stored value.
        wvf.nSpatialSamples = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'refpupilplanesize', 'refpupilplanesizemm', 'fieldsizemm'}
        % Total size of computed field in pupil plane.  This is for the measurement
        % wavelength.  The value can vary with wavelength, but this one
        % sets the scale for all the other wavelengths.
        %
        % This is a stored value.
        wvf.refSizeOfFieldMM = val;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'refpupilplanesampleinterval', 'refpupilplanesampleintervalmm', 'fieldsamplesize','fieldsamplesizemmperpixel'}
        % Pixel sampling interval of sample pupil field.  This is for the measurement
        % wavelength.  The value can vary with wavelength, but this one
        % sets the scale for all the other wavelengths.
        wvf.refSizeOfFieldMM = val*wvf.nSpatialSamples;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    case {'refpsfarcminpersample', 'refpsfarcminperpixel'}
        % Arc minutes per pixel of the sampled psf at the measurement
        % wavelength.
        %
        % When we convert between the pupil function and the PSF,
        % we use the fft.  Thus the size of the image in pixels
        % is the same for the sampled pupil function and the sampled
        % psf.
        %
        % The number of arc minutes per pixel in the sampled PSF is
        % related to the number of mm per pixel for hte pupil function,
        % with the relation depending on the wavelength.  The fundamental
        % formula in the pupil plane is that the pixel sampling interval
        % in cycles/radian is:
        %
        %   pupilPlaneCyclesRadianPerPix = pupilPlaneField/[lambda*npixels]
        %
        % where npixels is the number of linear pixels and lambda is the
        % wavelength. This formula may be found as Eq 10 of Ravikumar et al.
        % (2008), "Calculation of retinal image quality for polychromatic light,"
        % JOSA A, 25, 2395-2407, at least if we think their quantity d is the
        % size of the pupil plane field being sampled.
        %
        % If we now remember how units convert when we do the fft, we obtain
        % that the number of radians in the PSF image is the inverse of the
        % sampling interval:
        %
        %   radiansInPsfImage = [lambda*npixels]/pupilPlaneField
        %
        % which then gives us the number of radiansPerPixel in the
        % PSF image as
        %
        %   radiansPerPixel = lambda/pupilPlaneField
        %
        % The formula below implements this, with a conversion
        % from radians to minutes with factor (180*60/3.1416)
        % and converts wavelength to mm from nm with factor (.001*.001)
        %
        % DHB, 5/22/12, based on earler comments that were here.  Someone
        % else might take a look at the paper referenced above and the logic
        % of this comment and check that it all seems right.  Did I think
        % through the fft unit conversion correctly?  And, there must be
        % a more fundamental reference than the paper above, and for which
        % one wouldn't have to guess quite as much about what is meant.
        radiansPerPixel = val/(180*60/3.1416);
        wvf.refSizeOfFieldMM = wvfGet('measured wl','mm')/radiansPerPixel;
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
end

% Spectral
switch parm
    case {'calcwavelengths','wavelengths','wavelength','wls','wave'}
        % Normally just a vector of wavelengths in nm
        % but allow SToWls case for PTB users.
        %
        % wvfSet(wvfP,'wave',400:10:700)  OR
        % wvfSet(wvfP,'wave',[400 10 31])
        %
        % Note that it isn't sufficient just to call SToWls
        % (which will pass a vector of evenly spaced wls
        % through, because we might want to allow unevenly
        % spaced wls.)
        if size(val,2) == 3 && size(val,1) == 1 % SToWls case
            % Row vector with 3 entries.
            wls = SToWls(val);
            wvf.wls = MakeItWls(wls);
        else  % A column vector case
            wvf.wls = val(:);
        end
        PUPILFUNCTION_STALE = true;
        DIDASET = true;
        
    %case {'infocuswavelength','nominalfocuswl'}
    %    % In focus wavelength nm.  Single value.
    %    wvf.nominalFocusWl = val;
    case 'weightspectrum'
        % Used for calculating defocus to white light with many spectral
        % terms.
        if (length(val) ~= length(wvf.wls))
            error('Weighting spectrum dimension must match number of wavelengths');
        end
        wvf.weightingSpectrum = val;
        DIDASET = true;

        % Pupil parameters
    case {'calculatedpupil','calculatedpupildiameter'}
        % Pupil diameter in mm - must be smaller than measurements
        if (val > wvf.measpupilMM)
            error('Pupil diamter used for calculation must be smaller than that used for measurements');
        end
        wvf.calcpupilMM = val;
        DIDASET = true;

    case {'pupilfunction','pupilfunc'}
        % wvfSet(wvf,'pupil function',pf) - pf is a cell array of pupil
        % functions, one for each wavelength
        %
        % wvfSet(wvf,'pupil function',pf,[idx]) - pf pupil is for
        % wave(idx).
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
        DIDASET = true;
              
    case 'defocusdiopters'
        % Does not look like defocus is ever stored in diopters in other
        % functions. Only for user to set defocus diopters as an
        % alternative to using 4th term of zernike coeffs.
        wvf.defocusDiopters = val;           % Defocus
        DIDASET = true;

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
        DIDASET = true;
        
        % Special cases
    case {'sceparams','stilescrawford'}
        % The structure of sce is defined in sceCreate
        wvf.sceParams = val;
        DIDASET = true;
        
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
        DIDASET = true;
       
        %        % Not sure why this is set.  It is derived
        %     case 'strehl'
        %        %   wvf.strehl = val;
        
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
end

%% Catch the case where we don't know about the requested parameter
switch (parm)
    otherwise
        if (~DIDASET)
            error('Unknown parameter %s\n',parm);
        end
end

if (PUPILFUNCTION_STALE)
    
end

if (PSF_STALE)
    
end

return
