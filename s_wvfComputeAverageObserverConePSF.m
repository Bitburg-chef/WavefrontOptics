% wvfComputeAverageObserverConePSF
%
% Compute the cone PSFs of an average observer.  Done by combining measurements
% from a dataset of Zernike coefficients for a number of observers.
%
% This is done with respect to a specified spectral weighting function,
% and focus optimzed for a specified weighting of the different cone classes.
%
% See also: ComputeConePSFTest, wvfComputeConePSF
%   ComputePSFTest, wvfComputePSF, wvfComputePupilFunction,
%   GetStilesCrawfordParamsParams, wvfGetDefocusFromWavelengthDifference
%
% 8/29/11  dhb  Wrote it.

%% Clear
clear; close all;

%% Load cone sensitivities, set weighting spectrum.
S = [400 5 61];
wls = SToWls(S);
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);
load spd_D65
weightingSpectrum = SplineSpd(S_D65,spd_D65,S);

% Specify datafile for Zernike coefficients
zernikeFile = 'sampleZernikeCoeffs.txt';
measpupilMM = 6;
theZernikeCoeffs = load(zernikeFile);

wvfParams0.measpupilMM = measpupilMM;
wvfParams0.calcpupilMM = 3;
wvfParams0.wls = wls;
wvfParams0.nominalFocusWl = 550;
wvfParams0.defocusDiopters = 0;
wvfParams0.sizeOfFieldPixels = 201;
wvfParams0.sizeOfFieldMM = 16.212;
wvfParams0.T_cones = T_cones;
wvfParams0.weightingSpectrum = weightingSpectrum;
whichRow = floor(wvfParams0.sizeOfFieldPixels/2) + 1;

plotLimit = 2;
DOSCE = 0;
if (DOSCE)
    wvfParams0.sceParams = GetStilesCrawfordParams(wls,'berendshot');
else
    wvfParams0.sceParams = GetStilesCrawfordParams(wls,'none');
end
CIRCULARLYAVERAGE = 1;
coneWeights = [1 1 0];
criterionFraction = 0.9;

% Read coefficients and optimze PSF for each observer
for i = 1:size(theZernikeCoeffs,2)
    wvfParams = wvfParams0;
    wvfParams.zcoeffs = theZernikeCoeffs(:,i);
    wvfParams = wvfComputeOptimizedConePSF(wvfParams);
    arcminperpixel(i) = wvfParams.arcminperpixel;
    defocusDiopters(i) = wvfParams.defocusDiopters;
    lpsfo(:,:,i) = CenterPSF(wvfParams.conepsf(:,:,1));
    mpsfo(:,:,i) = CenterPSF(wvfParams.conepsf(:,:,2));
    spsfo(:,:,i) = CenterPSF(wvfParams.conepsf(:,:,3));
end

% Get optimized diffrac limited PSF
wvfParams = wvfParams0;
wvfParams.zcoeffs = zeros(61,1);
wvfParams = wvfComputeOptimizedConePSF(wvfParams);
arcminperpixeld = wvfParams.arcminperpixel;
defocusDioptersd = wvfParams.defocusDiopters;
lpsfd = CenterPSF(wvfParams.conepsf(:,:,1));
mpsfd = CenterPSF(wvfParams.conepsf(:,:,2));
spsfd = CenterPSF(wvfParams.conepsf(:,:,3));

% Get average LMS PSFs
avglpsfo = AverageOpticalPSFs(lpsfo);
avgmpsfo = AverageOpticalPSFs(lpsfo);
avgspsfo = AverageOpticalPSFs(lpsfo);
if (CIRCULARLYAVERAGE)
    avglpsfo = CircularlyAveragePSF(avglpsfo);
    avgmpsfo = CircularlyAveragePSF(avgmpsfo);
    avgspsfo = CircularlyAveragePSF(avgspsfo);
    lpsfd = CircularlyAveragePSF(lpsfd);
    mpsfd = CircularlyAveragePSF(mpsfd);
    spsfd = CircularlyAveragePSF(spsfd);
end
onedLPSFo = avglpsfo(whichRow,:);
onedMPSFo = avgmpsfo(whichRow,:);
onedSPSFo = avgspsfo(whichRow,:);
onedLPSFd = lpsfd(whichRow,:);
onedMPSFd = mpsfd(whichRow,:);
onedSPSFd = spsfd(whichRow,:);

arcminutes = arcminperpixel(1)*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
figure; clf;
subplot(1,3,1); hold on
plot(arcminutes(index),onedLPSFo(index),'r','LineWidth',4);
plot(arcminutes(index),onedLPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
if (CIRCULARLYAVERAGE)
    title('Circularized L cone PSF');
else
    title('L cone PSF');
end
subplot(1,3,2); hold on
plot(arcminutes(index),onedMPSFo(index),'g','LineWidth',4);
plot(arcminutes(index),onedMPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
if (CIRCULARLYAVERAGE)
    title('Circularized M cone PSF');
else
    title('M cone PSF');
end
subplot(1,3,3); hold on
plot(arcminutes(index),onedSPSFo(index),'b','LineWidth',4);
plot(arcminutes(index),onedSPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
if (CIRCULARLYAVERAGE)
    title('Circularized S cone PSF');
else
    title('S cone PSF');
end
drawnow;

% Save
if (CIRCULARLYAVERAGE)
    outfile = sprintf('AverageConePSFCIRC_%d_%d_%d_%0.1f',coneWeights(1),coneWeights(2),coneWeights(3),criterionFraction);
else
    outfile = sprintf('AverageConePSF_%d_%d_%d_%0.1f',coneWeights(1),coneWeights(2),coneWeights(3),criterionFraction);
end
save(outfile,'avglpsfo','avgmpsfo','avgspsfo','lpsfd','mpsfd','spsfd','arcminperpixel','defocusDiopters');





