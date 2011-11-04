function wvfP = wvfCreate
% Create the default wavefront parameters
%
%
%
%

% Wavelength samples
S = [400 5 61];
wls = SToWls(S);
% Something about the cones.  
T = load('T_cones_ss2');   % Probably in the PTB
T_cones = SplineCmf(T.S_cones_ss2,T.T_cones_ss2,S);
% vcNewGraphWin; plot(wave,T_cones'); xlabel('Wavelength (nm)');

% Weighting spectrum, probably for combining the PSFs in an average 
T = load('spd_D65');
weightingSpectrum = SplineSpd(T.S_D65,T.spd_D65,S);

% Parameters in wvf
wvfP.zcoeffs     = zeros(65,1);     % Zernicke coefficients
wvfP.measpupilMM = 6;               % Default pupil diameter?
wvfP.calcpupilMM = 3;               % Default pupil??? radius?
wvfP.wls         = wls;             % Wavelength samples
wvfP.nominalFocusWl = 550;          % In focus wavelength (nm)
wvfP.defocusDiopters = 0;           % Defocus
wvfP.sizeOfFieldPixels = 201;       % In pixels?  No units?
wvfP.sizeOfFieldMM = 16.212;        % Not sure which field
wvfP.T_cones = T_cones;             % Resampled cone spectral absorptions
wvfP.weightingSpectrum = weightingSpectrum;  % Probably used for combined psf
wvfP.sceParams = sceGetParams(wls,'berendshot');  % See notes therein
% This is an alternative:
% wvfP.sceParams = sceGetParams(wls,'none');

return
