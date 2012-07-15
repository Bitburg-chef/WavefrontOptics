% v_wvfComputeConePSF
%
% The the routines that compute L, M, and S cone PSFs from Zernike
% coefficients.
%
% See also: wvfComputeConePSF, wvfComputePSF, wvfComputePupilFunction,
%   sceGetParams, wvfGetDefocusFromWavelengthDifference
%
% The circular averaging is not a good idea for a single subject, but if
% you want to obtain an average over subjects it seems like a good idea.
%
% 8/21/11  dhb  Wrote it.
% 3/15/12  mdl  Edited to use wvfSet. Also updated to use fieldSampleSizeMMperPixel

%% Initialize
clear; close all;
s = which('v_wvfComputeConePSF');
cd(fileparts(s));
s_initISET;

%% Parameters
whichSubject = 1;
DOSCE = 0;
plotLimit = 2;
CIRCULARLYAVERAGE = 1;

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
theZernikeCoeffs = load('sampleZernikeCoeffs.txt','-ascii');
if (size(theZernikeCoeffs,2) == 585)
    theZernikeCoeffs = reshape(theZernikeCoeffs,9,65)';
end
if (size(theZernikeCoeffs,1) ~= 65 || size(theZernikeCoeffs,2) ~= 9)
    error('Surprising size for read in Hofer zSamples.')
end
wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'zcoeffs',theZernikeCoeffs(:,whichSubject));
wvf0 = wvfSet(wvf0,'calc wavelengths',wls);
if (DOSCE == 1)
    sce = sceCreate(wls,'berendshot');
    wvf0 = wvfSet(wvf0,'sce params',sce);
else
    sce = sceCreate(wls,'none');
    wvf0 = wvfSet(wvf0,'sce params',sce);
end

% Compute LMS psfs both for a subject and diffraction limited
wvfParams1 = wvf0;
wvfParams1 = wvfComputeConePSF(wvfParams1);
wvfParams2 = wvf0;
wvfParams2.zcoeffs = diffracZcoeffs;
wvfParams2 = wvfComputeConePSF(wvfParams2);

lpsf = psfCenter(wvfParams1.conepsf(:,:,1));
mpsf = psfCenter(wvfParams1.conepsf(:,:,2));
spsf = psfCenter(wvfParams1.conepsf(:,:,3));
lpsfd = psfCenter(wvfParams2.conepsf(:,:,1));
mpsfd = psfCenter(wvfParams2.conepsf(:,:,2));
spsfd = psfCenter(wvfParams2.conepsf(:,:,3));

if (CIRCULARLYAVERAGE)
    lpsf = psfCircularlyAverage(lpsf);
    mpsf = psfCircularlyAverage(mpsf);
    spsf = psfCircularlyAverage(spsf);
    lpsfd = psfCircularlyAverage(lpsfd);
    mpsfd = psfCircularlyAverage(mpsfd);
    spsfd = psfCircularlyAverage(spsfd);
end
whichRow = floor(wvfGet(wvfParams1,'npixels')/2) + 1;
arcminutes = wvfParams1.arcminperpix*((1:wvfGet(wvfParams1,'npixels'))-whichRow);

% Make a plot through the peak of the returned PSFs.
theFig = figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
set(gcf,'Position',position);
subplot(1,3,1); hold on
onedLPSF = lpsf(whichRow,:);
onedLPSFD = lpsfd(whichRow,:);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedLPSF(index),'r','LineWidth',2);
plot(arcminutes(index),onedLPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');

if (CIRCULARLYAVERAGE)
    title('Circularized L cone PSF');
else
    title('L cone PSF');
end

subplot(1,3,2); hold on
onedMPSF = mpsf(whichRow,:);
onedMPSFD = mpsfd(whichRow,:);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedMPSF(index),'g','LineWidth',2);
plot(arcminutes(index),onedMPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');

if (CIRCULARLYAVERAGE)
    title('Circularized M cone PSF');
else
    title('M cone PSF');
end

subplot(1,3,3); hold on
onedSPSF = spsf(whichRow,:);
onedSPSFD = spsfd(whichRow,:);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedSPSF(index),'b','LineWidth',2);
plot(arcminutes(index),onedSPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');

if (CIRCULARLYAVERAGE)
    title('Circularized S cone PSF');
else
    title('S cone PSF');
end
drawnow;

%% TEST2.  Optimize focus and add to the plot.
%
% This takes a long time.
wvfParams3 = wvf0;
wvfParams3.coneWeights = [1 1 0];
wvfParams3.criterionFraction = 0.9;
wvfParams3 = wvfComputeOptimizedConePSF(wvfParams3);
lpsfo = psfCenter(wvfParams3.conepsf(:,:,1));
mpsfo = psfCenter(wvfParams3.conepsf(:,:,2));
spsfo = psfCenter(wvfParams3.conepsf(:,:,3));
if (CIRCULARLYAVERAGE)
    lpsfo = psfCircularlyAverage(lpsfo);
    mpsfo = psfCircularlyAverage(mpsfo);
    spsfo = psfCircularlyAverage(spsfo);
end
onedLPSFo = lpsfo(whichRow,:);
onedMPSFo = mpsfo(whichRow,:);
onedSPSFo = spsfo(whichRow,:);

figure(theFig);
subplot(1,3,1);
plot(arcminutes(index),onedLPSFo(index),'r','LineWidth',4);
subplot(1,3,2); hold on
plot(arcminutes(index),onedMPSFo(index),'g','LineWidth',4);
subplot(1,3,3); hold on
plot(arcminutes(index),onedSPSFo(index),'b','LineWidth',4);
drawnow;





