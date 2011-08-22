% ComputePSFFromZernikeBasicTest
%
% Performs basic test of the routines that compute polychromatic PSFs from
% Zernike coefficients.
%
% 8/21/11  dhb  Wrote it, based on code provided by Heidi Hofer.

%% Clear
clear; close all;

%% TEST1: The Zernike code should return the diffraction limited PSF if we pass
% all zeros as the coefficients.  So the first test is whether this works,
% when we compare to what we get when we produce the diffraction limited
% PSF using analytic formulae implemented in the PTF routine AiryPattern.
%
% This appears to work correctly.
theWavelength = 650;
wavelengthOffset = 250;
pupilOffset = 4;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
diffracZcoeffs = zeros(65,1);
measpupilMM = 8;
calcpupilMM = 3;
defocusDiopters = 0;

% Make a plot through the peak of the returned PSF, normalized to peak of 1.
% Compare to what we get from PTB AiryPattern function -- should match
figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
set(gcf,'Position',position);
subplot(1,3,1); hold on
[diffracPSF1,arcminperpix] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,theWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
whichRow = floor(sizeOfFieldPixels/2) + 1;
onedPSF1 = diffracPSF1(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,calcpupilMM,theWavelength);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',calcpupilMM,theWavelength));

subplot(1,3,2); hold on
[diffracPSF1,arcminperpix] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength-wavelengthOffset,theWavelength-wavelengthOffset,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
whichRow = floor(sizeOfFieldPixels/2) + 1;
onedPSF1 = diffracPSF1(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,calcpupilMM,theWavelength-wavelengthOffset);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',calcpupilMM,theWavelength-wavelengthOffset));

subplot(1,3,3); hold on
[diffracPSF1,arcminperpix] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM+pupilOffset,theWavelength,theWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
whichRow = floor(sizeOfFieldPixels/2) + 1;
onedPSF1 = diffracPSF1(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,calcpupilMM+pupilOffset,theWavelength);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',calcpupilMM+pupilOffset,theWavelength));

%% TEST2: If we change the nominal focus away from a specified wavelength, the psf should
% get broader.  How much broader, I don't know but we can at least verify the qualitative
% behavior, again for the diffaction limited case.
%
% Note from DHB.  The psfs do get broader, and they take on multiple peaks.  Is this right?  I don't know.
nominalFocusWavelength = 550;
theWavelength = 550;
wavelengthOffset = 50;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 3;
defocusDiopters = 0;

figure; clf; hold on
[diffracPSF1,arcminperpix] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
whichRow = floor(sizeOfFieldPixels/2) + 1;
onedPSF1 = diffracPSF1(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
[diffracPSF2] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength-wavelengthOffset,nominalFocusWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
onedPSF2 = diffracPSF2(whichRow,:);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
[diffracPSF3] = ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength+wavelengthOffset,nominalFocusWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM);
onedPSF3 = diffracPSF3(whichRow,:);
plot(arcminutes(index),onedPSF3(index),'g','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title(sprintf('Effect of defocus, +/- %0.1f nm',wavelengthOffset));

%% TEST3.  Include the Stiles-Crawford effect.
%
% Note from DHB.  I'm not sure here what the effect should look like,
% so this mostly just demonstrates that the code runs.
nominalFocusWavelength = 550;
theWavelength = 550;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 8;
defocusDiopters = 0;
sceParams = GetStilesCrawfordParams(theWavelength,'berendshot');

figure; clf; hold on
[diffracPSF1,arcminperpix] = ...
    ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,[]);
whichRow = floor(sizeOfFieldPixels/2) + 1;
onedPSF1 = diffracPSF1(whichRow,:);
arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 4);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
[diffracPSF2,nil,strehlSCE,sceFrac] = ...
    ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
diffracPSF2 = CenterPSF(diffracPSF2);
onedPSF2 = diffracPSF2(whichRow,:);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title(sprintf('SCE, returned strehl %0.3f, direct from non-SCE diffrac %0.3f, frac. abs. is %0.2f',strehlSCE,max(diffracPSF2(:))/max(diffracPSF1(:)),sceFrac));

%% TEST4.  We have a set of measured Zernike coefficients for 9 subjects.  Let's make sure
% we can compute the psfs using these, and look at a slice through each of them.
%
% The computed PSFs are recentered, with their maximum in the center, so that we
% see the real peak when we take the 1D slice.
%
% This also optimizes the defocus to maximize the strehl ratio for each subject.
%
% Note from DHB.  Again, I don't know if these are correct, but at least you can see
% that you get a wide range of answers by using different subjects' data.
nominalFocusWavelength = 550;
theWavelength = 550;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 6;
calcpupilMM = 6;
defocusDiopters = 0;
theZernikeCoeffs = load('sampleZernikeCoeffs.txt');
DOSCE = 1;
if (DOSCE)
    sceParams = GetStilesCrawfordParams(theWavelength,'berendshot');
else
    sceParams = GetStilesCrawfordParams(theWavelength,'none');
end

figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
position(4) = 1600;
set(gcf,'Position',position);
for i = 1:9
    subplot(3,3,i); hold on
    [diffracPSF1,arcminperpix] = ...
        ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,[]);
    [diffracPSF2] = ...
        ComputePSFFromZernike(diffracZcoeffs,measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
    whichRow = floor(sizeOfFieldPixels/2) + 1;
    onedPSF1 = diffracPSF1(whichRow,:);
    arcminutes = arcminperpix*((1:sizeOfFieldPixels)-whichRow);
    index = find(abs(arcminutes) < 6);
    plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
    [thePSF2,nil,strehl2] = ...
        ComputePSFFromZernike(theZernikeCoeffs(:,i),measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,0,sizeOfFieldPixels,sizeOfFieldMM,sceParams);
    thePSF2 = CenterPSF(thePSF2);
    onedPSF2 = thePSF2(whichRow,:);
    plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
    strehlDirect2 = max(thePSF2(:))/max(diffracPSF1(:));
    strehlDirect2SCE = max(thePSF2(:))/max(diffracPSF2(:));
    
    % Optimize strehl by varying defocus
    bestStrehl = 0;
    defocusDiopters = -2:0.25:2;
    for j = 1:length(defocusDiopters)
        [thePSF3,nil,strehl3] = ...
        ComputePSFFromZernike(theZernikeCoeffs(:,i),measpupilMM,calcpupilMM,theWavelength,nominalFocusWavelength,defocusDiopters(j),sizeOfFieldPixels,sizeOfFieldMM,sceParams);
        if (strehl3 > bestStrehl)
            bestStrehl = strehl3;
            bestPSF3 = thePSF3;
            bestDefocusDiopters = defocusDiopters(j);
        end
    end
    bestPSF3 = CenterPSF(bestPSF3);
    onedPSF3 = bestPSF3(whichRow,:);
    plot(arcminutes(index),onedPSF3(index),'g','LineWidth',2);
    strehlDirect3 = max(bestPSF3(:))/max(diffracPSF1(:));
    strehlDirect3SCE = max(bestPSF3(:))/max(diffracPSF2(:));

      
    xlabel('Arc Minutes');
    ylabel('PSF');
    title(sprintf('Subject %d, strehl %0.2f (no defocus), %0.2f (defocus of %0.2f D)',i,strehl2,bestStrehl,bestDefocusDiopters));
    drawnow;
    
    if (DOSCE)
        fprintf('Subject %i, with SCE correction\n',i);
    else
        fprintf('Subject %i, no SCE correction\n',i);
    end
    fprintf('\tNo defocus: direct strehl with diffrac %0.3f, direct strehl with SCE diffrac %0.2f, returned, %0.2f\n',strehlDirect2,strehlDirect2SCE,strehl2);
    fprintf('\twith defocus: direct strehl with diffrac %0.3f, direct strehl with SCE diffrac %0.2f, returned, %0.2f\n',strehlDirect3,strehlDirect3SCE,bestStrehl);
        
end





