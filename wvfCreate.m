function wvfP = wvfCreate(varargin)
% Create the wavefront parameters structure
%
%    wvfP = wvfCreate(varargin)
%
% varargin:  Structured as param, val pairs
%   
% See also: wvfSet, wvfGet, sceCreate, sceGet
%
% Examples:
%    wvfP = wvfCreate('wave',[400:10:700]);
%
% (c) Wavefront Toolbox Team 2011, 2012

%% Book-keeping
wvfP.name = 'default';
wvfP.type = 'wvf';

%% Zernike Coefficients
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
wvfP.zcoeffs = zeros(65,1);      % Zernike coefficients
wvfP.measpupilMM = 8;            % Pupil diameter for measurements (mm)
wvfP.measWlNM = 550;             % Measurement wavelength (nm)
wvfP.measOpticalAxisDeg = 0;     % Measurement optical axis, degrees eccentric from fovea.
wvfP.measObserverAcommodationDiopters = 0; % Observer accommodation, in diopters relative to relaxed state of eye;

%% Fundamental sampling parameters
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
wvfP.constantSampleIntervalDomain = 'psf';  % Options are {'psf','pupil'}
wvfP.nSpatialSamples = 201;                 % Number of linear spatial samples
wvfP.pupilPlaneReferenceSizeMM = 16.212;    % Size of sampled pupil plane, referred to measurement wavelength

% Keep old code from breaking for now.  These should go away.
wvfP.sizeOfFieldMM = wvfP.refPupilPlaneSizeMM;
wvfP.fieldSampleSizeMMperPixel = wvfP.pupilPlaneReferenceSizeMM/wvfP.nSpatialSamples;

%% What to calculate for

% We can calculate the pupil function for any pupil diameter smaller
% than the diameter over which the measurements extend.  This defines
% the size to be used for the calculations represented by the wvfP
% object.
wvfP.calcpupilMM = 3;               % Used for this calculation

% Vector of wavelength samples to calclate for - Default is 550, monochromatic
wvfP.calcwlsNM = 550; 
wvfP.wls = wvfP.calcWlsNM;  % This is the old name, should go away eventually.      

% These are used in order to adjust the defocus term of the Zernike coeffs
% (zcoeffs(4)). Nominal Focus wavelength is compared to wavelength
% specified in the PSF calculation to find additional defocus.
% Used in wvfGetDefocusFromWavelengthDifference
wvfP.nominalFocusWl = 550;          % In focus wavelength (nm)
wvfP.defocusDiopters = 0;           % Defocus

% Something about the cones.  
% S is a length 3 vector of the format: [start spacing Nsamples]
% ex: S = [400 50 5]; 5 wavelength samples 400, 450, 500, 550, 600
S = [550 1 1]; 
T = load('T_cones_ss2');   % Probably in the PTB
T_cones = SplineCmf(T.S_cones_ss2,T.T_cones_ss2,S);
% vcNewGraphWin; plot(wave,T_cones'); xlabel('Wavelength (nm)');

% Weighting spectrum, for combining the PSFs in an average 
T = load('spd_D65');
weightingSpectrum = SplineSpd(T.S_D65,T.spd_D65,S);

wvfP.T_cones = T_cones;                      % Resampled cone spectral absorptions
wvfP.weightingSpectrum = weightingSpectrum;  % Probably used for combined psf

% Sets up the Stiles Crawford Effect parameters. 
wvfP.sceParams = sceCreate([],'none');

% Handle any additional arguments via wvfSet
if ~isempty(varargin)
    if isodd(length(varargin))
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        wvfP = wvfSet(wvfP,varargin{ii},varargin{ii+1});
    end
end

return
