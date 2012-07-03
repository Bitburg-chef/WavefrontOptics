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
%
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
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'2dpsf space','um',waveIdx,maxUM);

% Make a graph of the PSF within 2 arc min
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'2dpsf angle','min',waveIdx,maxMIN);

%% Plot the middle row of the psf, scaled to peak of 1
vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN);
hold on

% Used for plotting comparisons below
arcminutes = wvfGet(wvfParams,'support arcmin','min',waveIdx);
arcminpersample = wvfGet(wvfParams,'ref psf sample interval');
arcminpersample1 = wvfGet(wvfParams,'psf arcmin per sample',1);
arcminpersample2 = wvfGet(wvfParams,'psf angle per sample',[],1);
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

%% Repeat the calculation with a wavelength offset.  To keep
% the new wavelength in focus in the calculations, we add an
% explicit observer focus correction, with the amount computed
% by wvfLCAFromWavelengthDifference relative to the measured wavelength
newWave = 400;  
wvfParams1 = wvfParams0;
wvfParams1 = wvfSet(wvfParams1,'wave',newWave);
lcaDiopters = wvfLCAFromWavelengthDifference(newWave,wvfGet(wvfParams1,'measured wl'));
wvfParams1 = wvfSet(wvfParams1,'calc observer focus correction',lcaDiopters);
wvfParams = wvfComputePSF(wvfParams1);

vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN);
 
hold on
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

% PSF angular sampling should be the same across wavelengths
arcminpersample2 = wvfGet(wvfParams,'psf angle per sample',[],1);
if (arcminpersample2 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end

%% Repeat the calculation with a different pupil size at original wavelength
pupilMM = 7; 
wvfParams2 = wvfParams0;
wvfParams2 = wvfSet(wvfParams2,'calc pupil size',pupilMM);
wvfParams = wvfComputePSF(wvfParams2);

vcNewGraphWin([],'upper left');
wvfPlot(wvfParams,'1d psf angle normalized','min',waveIdx,maxMIN);
hold on;

onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));



