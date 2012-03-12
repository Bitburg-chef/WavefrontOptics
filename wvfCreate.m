function wvfP = wvfCreate(varargin)
% Create the wavefront parameters structure
%
%    wvfP = wvfCreate(varargin)
%
% varargin:  Will be structured as param, val pairs
%   
% See also:  scePGet
%
% Examples:
%    wvfP = wvfCreate('wave',[400:10:700]);
%    wvfP = wvfCreate('
%
% (c) WVF Toolbox Team 2011

% Should have a switch for reading varargin param, val pairs.
% Param reading/setting seems to be handled in wvfSet instead (KP 3/11/12)

% Wavelength samples - Default is 550, monochromatic
S = [550 1 1]; 
% Considering removing this and specifying only 1 wl at a time (KP)

% Book-keeping
wvfP.name = 'default';
wvfP.type = 'wvf';

% Parameters in wvf
wvfP.zcoeffs     = zeros(65,1);     % Zernike coefficients
%Zernike coeffs 0,1,2 (piston, verticle tilt, horizontal tilt) are
%typically 0 since they are either constant (piston) or only change the 
%point spread location, not quality, as measured in wavefront aberrations.
%65 coefficients represent the "j" single-index scheme of OSA standards
%for 10 orders of terms (each order has order+1 number of terms).
%10 orders, including 0th order, should result in 66 total terms, but 0th
%order term, coeff(0)/piston, is neglected, resulting in 65 terms.
%Thus, sample data typically starts with coeff(3), which is +/-45 degree
%astigmatism. (KP)

% These are diameters. The measpupilMM is the pupil diameter when the
% aberrations were measured.  The reason for specifying the diameter at
% measurement time is to allow an error check that you're not trying to
% compute for a pupil larger than that for which you have data. (DHB)
wvfP.measpupilMM = 8;               % Maximum pupil diameter
wvfP.calcpupilMM = 3;               % Used for this calculation

% These are used in order to adjust the defocus term of the Zernike coeffs
% (zcoeffs(4)). Nominal Focus wavelength is compared to wavelength
% specified in the PSF calculation to find additional defocus.
% Used in wvfGetDefocusFromWavelengthDifference
wvfP.nominalFocusWl = 550;          % In focus wavelength (nm)
wvfP.defocusDiopters = 0;           % Defocus

% What field are we talking about here?
% KP: Square field over which pupil function is computed, typically
% larger than the calculated pupil MM
wvfP.sizeOfFieldPixels = 201;       % In pixels?  No units?
wvfP.sizeOfFieldMM = 16.212;        % Not sure which field

% Should this be just 1 wavelength or a vector? DHB seems to use wls as the
% wavelength list and S as a 3-vector of [start sample number]

% Specifies the wavelength(s) over which the pupil function and PSF will 
% be calculated. Needed in order to find axial aberration (LCA) caused
% by not necessarily using nominal focus wavelength.
% Will modify zcoeffs(4) "Defocus" term.
wls = SToWls(S);
wvfP.wls         = wls;             % Wavelength samples

% Something about the cones.  
T = load('T_cones_ss2');   % Probably in the PTB
T_cones = SplineCmf(T.S_cones_ss2,T.T_cones_ss2,S);
% vcNewGraphWin; plot(wave,T_cones'); xlabel('Wavelength (nm)');

% Weighting spectrum, for combining the PSFs in an average 
T = load('spd_D65');
weightingSpectrum = SplineSpd(T.S_D65,T.spd_D65,S);

wvfP.T_cones = T_cones;             % Resampled cone spectral absorptions
wvfP.weightingSpectrum = weightingSpectrum;  % Probably used for combined psf

% Sets up the Stiles Crawford Effect parameters. 
wvfP.sceParams = sceCreate([],'none');

% Set additional arguments
if ~isempty(varargin)
    if isodd(length(varargin))
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        wvfP = wvfSet(wvfP,varargin{ii},varargin{ii+1});
    end
end

return
