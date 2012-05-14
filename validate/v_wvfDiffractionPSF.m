% v_wvfDiffractionPSF
%
% Tests the monochromatic PSFs computed from Zernike coefficients. The
% script compares the computations with those in PTB.
%
% See also: wvfComputePSF, wvfComputePupilFunction,
%           sceCreate, wvfGetDefocusFromWavelengthDifference
%
% TODO
%   a) Add wvfPlot to output normalized to 1 curve.
%   b) Compare with ISET.  The implementation there uses spatial units, not
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
% (c) Wavefront Toolbox Team, 2012

%% Clear
% clear; close all;

%% Or
s_initISET

%% Compare pointspread function in wvf with psf in Psych Toolbox

% When the Zernike coefficients are all zero, the wvfComputePSF code should
% return the diffraction limited PSF.  We test whether this works by
% comparing to the diffraction limited PSF implemented in the PTB routine
% AiryPattern.

% Set up parameters structure
wvfParams0 = wvfCreate;

% Plotting ranges for MM, UM, and Minutes of angle
maxMM = 1;
maxUM = 20;
maxMIN = 2;

% Which wavelength index (wave(idx)) (there is only one) to plot
waveIdx = 1;

% Calculate the PSF, normalized to peak of 1.
wvfParams = wvfComputePSF(wvfParams0);

% Make a graph of the PSF within 1 mm of center
vcNewGraphWin;
wvfPlot(wvfParams,'2dpsf space','um',waveIdx,maxUM);

% Make a graph of the PSF within 2 arc min
vcNewGraphWin;
wvfPlot(wvfParams,'2dpsf angle','min',waveIdx,maxMIN);

%% Plot the middle row of the psf, scaled to peak of 1
vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN);
hold on

% Used for plotting comparisons below
arcminutes = wvfGet(wvfParams,'support arcmin');
index = find(abs(arcminutes) < 2);
radians = (pi/180)*(arcminutes/60);

% Compare to what we get from PTB AiryPattern function -- should match
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% Repeat the calculation with a wavelength offset

% OK, we seem to have a problem.  Probably because I didn't do the right
% correction in wvfComputePupilFunction.  See the notes in there.  The
% stuff only runs right at 550, not at other wavelengths.

newWave = 400;  
wvfParams1 = wvfParams0;
wvfParams1 = wvfSet(wvfParams1,'wave',newWave);
wvfParams1 = wvfSet(wvfParams1,'in focus wavelength', newWave);
wvfParams = wvfComputePSF(wvfParams1);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN)
 
hold on
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% Repeat the calculation with a different pupil size
pupilMM = 7;   % In millimeters?
wvfParams2 = wvfParams0;
wvfParams2 = wvfSet(wvfParams2,'calculated pupil',pupilMM);
wvfParams = wvfComputePSF(wvfParams2);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN);
hold on;

onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% Put ISET diffraction comparisons here or in the next script


