%% v_wvfStilesCrawford
%
% Note from DHB.  I'm not sure here what the effect should look like, given
% that we're using a diffraction limited pupil, so this mostly just
% demonstrates that the code runs.
%
% Note from HH: This looks about right.

%%
wvfParams0 = wvfCreate;

% In focus wavelength and other parameters
nominalFocusWavelength = 550;
theWavelength = 550;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 8;  % Bigger pupil to see the effect
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

wvfParams0 = wvfSet(wvfParams0,'sce params',sceCreate(theWavelength,'berendshot'));

%% No Stiles Crawford effect

wvfParams = wvfSet(wvfParams0,'sce params',[]);
wvfParams = wvfComputePSF(wvfParams);

vcNewGraphWin; hold on
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);

% Make a graph of the PSF within 1 mm of center
% vcNewGraphWin;
% wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

%%  Leave the SCE in place

wvfParams = wvfComputePSF(wvfParams0);

% vcNewGraphWin;
% wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

wvfPlot(wvfParams,'1d psf angle','min',maxMIN);
hold on

strehl = wvfGet(wvfParams,'strehl');
title(sprintf('Strehl ratio:  %.1f',strehl));

% Figure this out
% title(sprintf('SCE, returned strehl %0.3f, direct from non-SCE diffrac %0.3f, frac. abs. is %0.2f',...
%     strehl,max(wvfParams2.psf(:))/max(wvfParams1.psf(:)),wvfParams2.sceFrac));

%% End
