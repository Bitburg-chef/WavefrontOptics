function wvf = wvfCreate(varargin)
% wvf = wvfCreate(varargin)
%
% Create the wavefront parameters structure
%
% varargin:  Structured as param, val pairs
%   
% See also: wvfSet, wvfGet, sceCreate, sceGet
%
% Examples:
%    wvf = wvfCreate('wave',[400:10:700]);
%
% (c) Wavefront Toolbox Team 2011, 2012

%% Book-keeping
wvf = [];
wvf = wvfSet(wvf,'name','default');
wvf = wvfSet(wvf,'type','wvf');

%% Zernike coefficients and related
wvf = wvfSet(wvf,'zcoeffs',zeros(65,1));
wvf = wvfSet(wvf,'measured pupil',8);
wvf = wvfSet(wvf,'measured wl',550);
wvf = wvfSet(wvf,'measured optical axis',0);
wvf = wvfSet(wvf,'measured observer accommodation',0);

%% Spatial sampling parameters
wvf = wvfSet(wvf,'sample interval domain','psf');
wvf = wvfSet(wvf,'spatial samples',201);
wvf = wvfSet(wvf,'ref pupil plane size',16.212);

%% Spectral
wvf = wvfSet(wvf,'calc wavelengths',550)

%% What to calculate for

% We can calculate the pupil function for any pupil diameter smaller
% than the diameter over which the measurements extend.  This defines
% the size to be used for the calculations represented by the wvf
% object.
wvf.calcpupilMM = 3;               % Used for this calculation
    

% These are used in order to adjust the defocus term of the Zernike coeffs
% (zcoeffs(4)). Nominal Focus wavelength is compared to wavelength
% specified in the PSF calculation to find additional defocus.
% Used in wvfGetDefocusFromWavelengthDifference
wvf.defocusDiopters = 0;           % Defocus

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

wvf.T_cones = T_cones;                      % Resampled cone spectral absorptions
wvf.weightingSpectrum = weightingSpectrum;  % Probably used for combined psf

% Sets up the Stiles Crawford Effect parameters. 
wvf.sceParams = sceCreate([],'none');

% Handle any additional arguments via wvfSet
if ~isempty(varargin)
    if isodd(length(varargin))
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        wvf = wvfSet(wvf,varargin{ii},varargin{ii+1});
    end
end

return
