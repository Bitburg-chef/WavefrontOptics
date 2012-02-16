% v_wvfDiffractionPSF
%
% Tests the monochromatic PSFs computed from Zernike coefficients. The
% script compares the computations with those in PTB.
%
% See also: wvfComputePSF, wvfComputePupilFunction,
%           sceCreate, wvfGetDefocusFromWavelengthDifference
%
% TODO
%   Compare with ISET.  The implementation there uses spatial units, not
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
clear; close all;

%% Compare with diffraction in PTB

% When the Zernike coefficients are all zero, the wvfComputePSF code should
% return the diffraction limited PSF.  We test whether this works by
% comparing to the diffraction limited PSF implemented in the PTB routine
% AiryPattern.

% Set up parameters structure
wvfParams0 = wvfCreate;

% Calculate the PSF, normalized to peak of 1.
wvfParams = wvfComputePSF(wvfParams0);

% Make a graph of the PSF within 1 mm of center
vcNewGraphWin;
maxMM = 1;
wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

% Make a graph of the PSF within 2 arc min
vcNewGraphWin;
maxMIN = 2;
wvfPlot(wvfParams,'2dpsf angle','min',maxMIN);

%% Plot the middle row of the psf, scaled to peak of 1
vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);
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

newWave = 400;  
wvfParams1 = wvfParams0;
wvfParams1 = wvfSet(wvfParams1,'wave',newWave);
wvfParams1 = wvfSet(wvfParams1,'in focus wavelength', newWave);
wvfParams = wvfComputePSF(wvfParams1);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN)
 
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
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);
hold on;

onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% Put ISET diffraction comparisons here or insert above

% ISET requires that we specify sizes.  So make a small scene and compute
% through to the optical image
scene = sceneCreate('lined65');
scene = sceneSet(scene,'hfov',0.2);
oi = oiCreate;
oi = oiCompute(scene,oi);

% For this scene we have fairly fine spatial angular resolution
% oiGet(oi,'angular resolution')*60  % In minutes
v = oiGet(oi,'angular support','min');
xAng = v(:,1,2);
yAng = v(1,:,1);

optics = oiGet(oi,'optics');
d = plotOTF(oi,'psf550');
