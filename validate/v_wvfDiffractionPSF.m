% v_wvfDiffractionPSF
%
% Tests the monochromatic PSFs computed from Zernike coefficients. The
% script compares the computations with those in PTB.
%
% See also: wvfCreate, wvfGet, wvfSet, wvfComputePSF, wvfComputePupilFunction,
%   wvfLCAFromWavelengthDifference
%
% TODO
%   a) Compare with ISET.  The implementation there uses spatial units, not
%   angles, to specify the samples. Further, it uses f-number, not just
%   pupil size.  So we need to build up ISET so that we can compute the psf
%   in terms of visual angle using only a pupil diameter (e.g., 3mm)
%   irrespective of the focal length (or equivalently f-number).
%
%   Is it the case that if we compute 
%     * the diffraction limited function for any f-number 
%     * use the focal length to  convert the spatial samples to arcmin 
%     * the result is the function for that aperture size
%
% 6/5/12  dhb  Got this to work with current code.
% 7/1/12  bw   Replaced waveIdx with wList in wvfGet and here and ...
%              Maybe we should just use GIT for these comments.
% (c) Wavefront Toolbox Team, 2012

%% Initialize
% cd(fileparts(mfilename('fullpath')));
s_initISET;

%% Compare pointspread function in wvf with psf in Psych Toolbox

% When the Zernike coefficients are all zero, the wvfComputePSF code should
% return the diffraction limited PSF.  We test whether this works by
% comparing to the diffraction limited PSF implemented in the PTB routine
% AiryPattern.

% Set up parameters structure
wvf0 = wvfCreate;

% Plotting ranges for MM, UM, and Minutes of angle
maxMM = 1;
maxUM = 20;
maxMIN = 2;

% Which wavelength index (wave(idx)) (there is only one) to plot
wList = wvfGet(wvf0,'wave');

% Calculate the PSF, normalized to peak of 1.
wvfParams = wvfComputePSF(wvf0);

% Make a graph of the PSF within maxUM of center
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'2dpsf space','um',wList,maxUM);

% Make a graph of the PSF within 2 arc min
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'2dpsf angle','min',wList,maxMIN);

%% Plot the middle row of the psf, scaled to peak of 1
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'1d psf angle normalized','min',wList,maxMIN);
hold on

% Used for plotting comparisons below
arcminutes = wvfGet(wvfParams,'support arcmin','min',wList);
arcminpersample = wvfGet(wvfParams,'ref psf sample interval');
arcminpersample1 = wvfGet(wvfParams,'psf arcmin per sample',wList);
arcminpersample2 = wvfGet(wvfParams,'psf angle per sample',[],wList);
if (arcminpersample1 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end
if (arcminpersample2 ~= arcminpersample1)
    error('Default units of get on ''psfanglepersample'' unexpectedly changed');
end

index = find(abs(arcminutes) < 2);
radians = (pi/180)*(arcminutes/60);

% Compare to what we get from PTB AiryPattern function -- should match
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% Repeat the calculation with a wavelength offset.  
% To keep the new wavelength in focus in the calculations, we add an
% explicit observer focus correction, with the amount computed by
% wvfLCAFromWavelengthDifference relative to the measured wavelength
wList = 400;  
wvf1 = wvf0;
wvf1 = wvfSet(wvf1,'wave',wList);
lcaDiopters = wvfLCAFromWavelengthDifference(wList,wvfGet(wvf1,'measured wl'));
wvf1 = wvfSet(wvf1,'calc observer focus correction',lcaDiopters);
wvf1 = wvfComputePSF(wvf1);

vcNewGraphWin([],'upper left');
w = wvfGet(wvf1,'wave');
pupilSize = wvfGet(wvf1,'calcpupilsize','mm');

wvfPlot(wvf1,'1d psf angle normalized','min',w,maxMIN);
hold on
onedPSF2 = AiryPattern(radians,pupilSize,w);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,w));

% PSF angular sampling should be the same across wavelengths
arcminpersample2 = wvfGet(wvf1,'psf angle per sample','min',w);
if (arcminpersample2 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end

%% Repeat the calculation with a different pupil size at original wavelength
pupilMM = 7; 
wList = 550;
wvf2  = wvf0;
wvf2  = wvfSet(wvf2,'wave',wList);
wvf2  = wvfSet(wvf2,'calc pupil size',pupilMM);

% lcaDiopters = wvfLCAFromWavelengthDifference(wList,wvfGet(wvf2,'measured wl'));
% wvf2 = wvfSet(wvf2,'calc observer focus correction',lcaDiopters);
wvf2.PSF_STALE = 1;
wvf2  = wvfComputePSF(wvf2);

wList = wvfGet(wvf2,'wave');
pupilSize = wvfGet(wvf2,'calcpupilsize','mm');

% Compare the curves
vcNewGraphWin([],'upper left');
wvfPlot(wvf2,'1d psf angle normalized','min',wList,maxMIN);
onedPSF2 = AiryPattern(radians,pupilSize,wList);

hold on
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,wList));

%%

