% v_wvfWaveDefocus
%
% Change the nominal focus away from a specified wavelength
%
% The psf should get broader. How much broader, I don't know but we can at
% least verify the qualitative behavior, again for the diffaction limited
% case.
%
% Note from DHB.  The psfs do get broader, and they take on multiple peaks.
% Is this right?  I don't know.
%
% Note from HH: The PSF should start to take on multiple peaks (ringing,
% etc should still be radially symmetric) with change in defocus, so this
% looks right.  At some point defocus will be large enough that you'll run
% into sampling problems.  You could check this by starting with very high
% sampling density and decreasing until just before the point where you
% start to see issues.
%
% See also:  v_wvfDiffractionPSF, v_wvfWaveDefocus
%
% (c) Wavefront Toolbox Team, 2012

%%

%s_initISET

%% Calculate point spread for wavelength defocus

% Set up default parameters structure with diffraction limited default
wvfParams0 = wvfCreate;

% In focus wavelength and other parameters
nominalFocusWavelength = 550;
theWavelength = 550;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 3;
defocusDiopters = 0;

% Set up parameters structure for basic wavefront
wvfParams0 = wvfSet(wvfParams0,'zcoeffs',zeros(65,1));
wvfParams0 = wvfSet(wvfParams0,'measured pupil',measpupilMM);
wvfParams0 = wvfSet(wvfParams0,'calculated pupil',calcpupilMM);
wvfParams0 = wvfSet(wvfParams0,'wave',theWavelength);
wvfParams0 = wvfSet(wvfParams0,'infocus wavelength',nominalFocusWavelength);
wvfParams0 = wvfSet(wvfParams0,'defocus diopters',defocusDiopters);

% This looks like it might be redundant and should be removed
wvfParams0 = wvfSet(wvfParams0,'field size pixels',sizeOfFieldPixels);
wvfParams0 = wvfSet(wvfParams0,'field size mm',sizeOfFieldMM);

% For plotting limits
maxMIN = 2;
maxMM  = 1;

%% Compute and plot the default
wvfParams = wvfComputePSF(wvfParams0);
vcNewGraphWin; hold on
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);

%% Create a wavelength offset of 50 nm for the best focus

% Offset the stimulus wavelength from the in focus wavelength
inFocus = wvfGet(wvfParams0,'inFocus Wavelength');
wavelengthOffset = 50;
wvfParams = wvfSet(wvfParams0,'wave',inFocus - wavelengthOffset);
wvfParams = wvfComputePSF(wvfParams);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);

vcNewGraphWin;
wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

%% Now shift the offset the other way

inFocus = wvfGet(wvfParams0,'inFocus Wavelength');
wavelengthOffset = 50;
wvfParams = wvfSet(wvfParams0,'wave',inFocus + wavelengthOffset);
wvfParams = wvfComputePSF(wvfParams);

wvfParams = wvfComputePSF(wvfParams);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);
title(sprintf('Effect of defocus, + %0.1f nm',wavelengthOffset));


vcNewGraphWin;
wvfPlot(wvfParams,'2dpsf space','mm',maxMM);


%% End