% ComputeLMSPSFFromZernikeTest
%
% Performs basic test of the routines that compute L, M, and S cone PSFs from
% Zernike coefficients.
%
% See also: ComputeLMSPSFFromZernike, ComputePSFFromZernike, ComputePupilFunctionFromZernike,
%   GetStilesCrawfordParams, GetDefocusFromWavelengthDifference
%
% The circular averaging is not a good idea for a single subject, but if you want
% to obtain an average over subjects it seems like a good idea.
%
% 8/21/11  dhb  Wrote it.

%% Clear
clear; close all;

%% Load cone sensitivities, set weighting spectrum.
S = [400 5 61];
wls = SToWls(S);
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);
load spd_D65
weightingSpectrum = SplineSpd(S_D65,spd_D65,S);

%% TEST1: Just make sure the code runs for some
% choice of parameters.  Focus is not optimized
% to produce the best PSFs.
%
% This appears to work correctly.
whichSubject = 1;
nominalFocusWl = 550;
defocusDiopters = 0;
wavelengthOffset = 250;
pupilOffset = 4;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
diffracZcoeffs = zeros(65,1);
theZernikeCoeffs = load('sampleZernikeCoeffs.txt');
zcoeffs = theZernikeCoeffs(:,whichSubject);
measpupilMM = 6;
calcpupilMM = 3;
plotLimit = 2;
DOSCE = 0;
if (DOSCE)
    sceParams = GetStilesCrawfordParams(wls,'berendshot');
else
    sceParams = GetStilesCrawfordParams(wls,'none');
end
CIRCULARLYAVERAGE = 1;

% Compute LMS psfs both for a subject and diffraction limited
[conepsf,arcminperpix] = ...
    ComputeLMSPSFFromZernike(wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
[conepsfd,arcminperpix] = ...
    ComputeLMSPSFFromZernike(wls,T_cones,weightingSpectrum,diffracZcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
lpsf = CenterPSF(conepsf(:,:,1));
mpsf = CenterPSF(conepsf(:,:,2));
spsf = CenterPSF(conepsf(:,:,3));
lpsfd = CenterPSF(conepsfd(:,:,1));
mpsfd = CenterPSF(conepsfd(:,:,2));
spsfd = CenterPSF(conepsfd(:,:,3));
if (CIRCULARLYAVERAGE)
    lpsf = CircularlyAveragePSF(lpsf);
    mpsf = CircularlyAveragePSF(mpsf);
    spsf = CircularlyAveragePSF(spsf); 
    lpsfd = CircularlyAveragePSF(lpsfd); 
    mpsfd = CircularlyAveragePSF(mpsfd); 
    spsfd = CircularlyAveragePSF(spsfd);
end
whichRow = floor(sizeOfFieldPixels/2) + 1;

% Make a plot through the peak of the returned PSFs.
figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
set(gcf,'Position',position);
subplot(1,3,1); hold on
onedLPSF = lpsf(whichRow,:);
onedLPSFD = lpsfd(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedLPSF(index),'r','LineWidth',2);
plot(arcminutes(index),onedLPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('L cone PSF');

subplot(1,3,2); hold on
onedMPSF = mpsf(whichRow,:);
onedMPSFD = mpsfd(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedMPSF(index),'g','LineWidth',2);
plot(arcminutes(index),onedMPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('M cone PSF');

subplot(1,3,3); hold on
onedSPSF = spsf(whichRow,:);
onedSPSFD = spsfd(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedSPSF(index),'b','LineWidth',2);
plot(arcminutes(index),onedSPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('S cone PSF');
drawnow;

%% TEST2.  Optimize focus and add to the plot.
%
% This takes a long time.
coneWeights = [1 1 0];
criterionFraction = 0.9;
[conepsfo,arcminperpixel,defocusDiopters] = ...
    ComputeOptimizedLMSPsfs(coneWeights,criterionFraction,wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,...
    sizeOfFieldPixels,sizeOfFieldMM,sceParams);
lpsfo = CenterPSF(conepsfo(:,:,1));
mpsfo = CenterPSF(conepsfo(:,:,2));
spsfo = CenterPSF(conepsfo(:,:,3));
if (CIRCULARLYAVERAGE)
    lpsfo = CircularlyAveragePSF(lpsfo);
    mpsfo = CircularlyAveragePSF(mpsfo);
    spsfo = CircularlyAveragePSF(spsfo); 
end
onedLPSFo = lpsfo(whichRow,:);
onedMPSFo = mpsfo(whichRow,:);
onedSPSFo = spsfo(whichRow,:);
subplot(1,3,1);
plot(arcminutes(index),onedLPSFo(index),'r','LineWidth',4);
subplot(1,3,2); hold on
plot(arcminutes(index),onedMPSFo(index),'g','LineWidth',4);
subplot(1,3,3); hold on
plot(arcminutes(index),onedSPSFo(index),'b','LineWidth',4);





