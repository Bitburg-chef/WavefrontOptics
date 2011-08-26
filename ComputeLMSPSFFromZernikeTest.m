% ComputeLMSPSFFromZernikeTest
%
% Performs basic test of the routines that compute L, M, and S cone PSFs from
% Zernike coefficients.
%
% See also: ComputeLMSPSFFromZernike, ComputePSFFromZernike, ComputePupilFunctionFromZernike,
%   GetStilesCrawfordParams, GetDefocusFromWavelengthDifference
%
% 8/21/11  dhb  Wrote it, based on code provided by Heidi Hofer.

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

% Compute LMS psfs both for a subject and diffraction limited
[lpsf,mpsf,spsf,arcminperpix] = ...
    ComputeLMSPSFFromZernike(wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
[lpsfd,mpsfd,spsfd,arcminperpix] = ...
    ComputeLMSPSFFromZernike(wls,T_cones,weightingSpectrum,diffracZcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
lpsf = CenterPSF(lpsf);
mpsf = CenterPSF(mpsf);
spsf = CenterPSF(spsf);
lpsfd = CenterPSF(lpsfd);
mpsfd = CenterPSF(mpsfd);
spsfd = CenterPSF(spsfd);
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
plot(arcminutes(index),onedLPSF(index),'r','LineWidth',4);
plot(arcminutes(index),onedLPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('L cone PSF');

subplot(1,3,2); hold on
onedMPSF = mpsf(whichRow,:);
onedMPSFD = mpsfd(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedMPSF(index),'g','LineWidth',4);
plot(arcminutes(index),onedMPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('M cone PSF');

subplot(1,3,3); hold on
onedSPSF = spsf(whichRow,:);
onedSPSFD = spsfd(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < plotLimit);
plot(arcminutes(index),onedSPSF(index),'b','LineWidth',4);
plot(arcminutes(index),onedSPSFD(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title('S cone PSF');
drawnow;

%% TEST4.  We have a set of measured Zernike coefficients for 9 subjects.  Let's make sure
% we can compute the psfs using these, and look at a slice through each of them.
%
% The computed PSFs are recentered, with their maximum in the center, so that we
% see the real peak when we take the 1D slice.
%
% This also optimizes the defocus to maximize the strehl ratio for each subject,
% so you can see the (large) effect of doing that.
%
% Note from DHB.  Again, I don't know if these are correct, but at least you can see
% that you get a wide range of answers by using different subjects' data.
% nominalFocusWavelength = 550;
% theWavelength = 550;
% sizeOfFieldPixels = 201;
% sizeOfFieldMM = 16.212;
% measpupilMM = 6;
% calcpupilMM = 3;
% defocusDiopters = 0;
% theZernikeCoeffs = load('sampleZernikeCoeffs.txt');
% DOSCE = 0;
% if (DOSCE)
%     sceParams = GetStilesCrawfordParams(theWavelength,'berendshot');
% else
%     sceParams = GetStilesCrawfordParams(theWavelength,'none');
% end
% 
% figure; clf;
% position = get(gcf,'Position');
% position(3) = 1600;
% position(4) = 1600;
% set(gcf,'Position',position);
% for i = 1:9
%     subplot(3,3,i); hold on
%     [diffracPSF1,arcminperpix] = ...
%         ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,[]);
%     [diffracPSF2] = ...
%         ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
%     whichRow = floor(sizeOfFieldPixels/2) + 1;
%     onedPSF1 = diffracPSF1(whichRow,:);
%     arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
%     index = find(abs(arcminutes) < 6);
%     plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
%     [thePSF2,nil,strehl2] = ...
%         ComputePSFFromZernike(theZernikeCoeffs(:,i),measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
%     thePSF2 = CenterPSF(thePSF2);
%     onedPSF2 = thePSF2(whichRow,:);
%     plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
%     strehlDirect2 = max(thePSF2(:))/max(diffracPSF1(:));
%     strehlDirect2SCE = max(thePSF2(:))/max(diffracPSF2(:));
%     
%     % Optimize strehl by varying defocus
%     bestStrehl = 0;
%     defocusDiopters = -2:0.25:2;
%     for j = 1:length(defocusDiopters)
%         [thePSF3,nil,strehl3] = ...
%         ComputePSFFromZernike(theZernikeCoeffs(:,i),measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,defocusDiopters(j),sizeOfFieldPixels,sizeOfFieldMM,sceParams);
%         if (strehl3 > bestStrehl)
%             bestStrehl = strehl3;
%             bestPSF3 = thePSF3;
%             bestDefocusDiopters = defocusDiopters(j);
%         end
%     end
%     bestPSF3 = CenterPSF(bestPSF3);
%     onedPSF3 = bestPSF3(whichRow,:);
%     plot(arcminutes(index),onedPSF3(index),'g','LineWidth',2);
%     strehlDirect3 = max(bestPSF3(:))/max(diffracPSF1(:));
%     strehlDirect3SCE = max(bestPSF3(:))/max(diffracPSF2(:));
% 
%     xlabel('Arc Minutes');
%     ylabel('PSF');
%     title(sprintf('Subject %d, strehl %0.2f (no defocus), %0.2f (defocus of %0.2f D)',i,strehl2,bestStrehl,bestDefocusDiopters));
%     drawnow;
%     
%     if (DOSCE)
%         fprintf('Subject %i, with SCE correction\n',i);
%     else
%         fprintf('Subject %i, no SCE correction\n',i);
%     end
%     fprintf('\tNo defocus: direct strehl with diffrac %0.3f, direct strehl with SCE diffrac %0.2f, returned, %0.2f\n',strehlDirect2,strehlDirect2SCE,strehl2);
%     fprintf('\twith defocus: direct strehl with diffrac %0.3f, direct strehl with SCE diffrac %0.2f, returned, %0.2f\n',strehlDirect3,strehlDirect3SCE,bestStrehl);
%         
% end





