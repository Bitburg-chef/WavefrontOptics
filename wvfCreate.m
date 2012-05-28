function wvfP = wvfCreate(varargin)
% wvfP = wvfCreate(varargin)
%
% Create the wavefront parameters structure
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

%% Zernike coefficients and related
wvfP.zcoeffs = zeros(65,1);      % Zernike coefficients
wvfP.measpupilMM = 8;            % Pupil diameter for measurements (mm)
wvfP.measWlNM = 550;             % Measurement wavelength (nm)
wvfP.measOpticalAxisDeg = 0;     % Measurement optical axis, degrees eccentric from fovea.
wvfP.measObserverAcommodationDiopters = 0; % Observer accommodation, in diopters relative to relaxed state of eye;

%% Spatial sampling parameters
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
