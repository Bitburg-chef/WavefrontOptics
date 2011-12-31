% s_wvfComputePSFTest
%
% Performs basic test of the routines that compute monochromatic PSFs from
% Zernike coefficients.
%
% See also: wvfComputePSF, wvfComputePupilFunction,
%           sceCreate, wvfGetDefocusFromWavelengthDifference
%
% 8/21/11  dhb  Wrote it, based on code provided by Heidi Hofer.
% 9/7/11   dhb  Got this working with wvfParams i/o.
%

% TODO
%   Consider plotting in terms of physical distance in the image plane,
%   rather than angle

% Include the WavefrontOpticsToolbox path and the Psychtoolbox path
%   addpath(genpath(pwd))
%   ptbPath

%% Clear
clear; close all;

%% Compare with diffraction
% The Zernike code should return the diffraction limited PSF if we pass all
% zeros as the coefficients.  So the first test is whether this works, when
% we compare to what we get when we produce the diffraction limited PSF
% using analytic formulae implemented in the PTB routine AiryPattern.
%
% This appears to work correctly.

% Set up parameters structure
wvfParams0 = wvfCreate;

% Calculate the PSF, normalized to peak of 1.
wvfParams = wvfComputePSF(wvfParams0);
% vcNewGraphWin; mesh(wvfParams.psf)

% Extract a row of the psf
whichRow = floor(wvfParams.sizeOfFieldPixels/2) + 1;
onedPSF1 = wvfParams.psf(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = wvfParams.arcminperpix*((1:wvfParams.sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);

% Make a plot through the peak of the returned PSF, normalized to peak of 1.
figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
set(gcf,'Position',position);
subplot(1,3,1); hold on
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);

% Compare to what we get from PTB AiryPattern function -- should match
radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

% Repeat the calculation with a wavelength offset
wavelengthOffset = 250;
wvfParams1 = wvfParams0;
wvfParams1.wls = wvfParams1.wls - wavelengthOffset;
wvfParams1.nominalFocusWl = wvfParams1.nominalFocusWl - wavelengthOffset;

[wvfParams] = wvfComputePSF(wvfParams1);
whichRow = floor(wvfParams.sizeOfFieldPixels/2) + 1;
onedPSF1 = wvfParams.psf(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = wvfParams.arcminperpix*((1:wvfParams.sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);

subplot(1,3,2); hold on
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));

plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

% Repeat the calculation with a pupil offset
pupilOffset = 4;   % In millimeters?
wvfParams2 = wvfParams0;
wvfParams2.calcpupilMM = wvfParams2.calcpupilMM + pupilOffset;
[wvfParams] = wvfComputePSF(wvfParams2);
whichRow = floor(wvfParams.sizeOfFieldPixels/2) + 1;
onedPSF1 = wvfParams.psf(whichRow,:);
onedPSF1 = onedPSF1/max(onedPSF1(:));
arcminutes = wvfParams.arcminperpix*((1:wvfParams.sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);

subplot(1,3,3); hold on
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);

radians = (pi/180)*(arcminutes/60);
onedPSF2 = AiryPattern(radians,wvfParams.calcpupilMM,wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfParams.calcpupilMM,wvfParams.wls(1)));

%% TEST2: Change the nominal focus away from a specified wavelength
%
% The psf should get broader. How much broader, I don't know but we can at
% least verify the qualitative behavior, again for the diffaction limited
% case.
%
% Note from DHB.  The psfs do get broader, and they take on multiple peaks.
% Is this right?  I don't know.
%
% Note from HH: The PSF should start to take on multiple peaks (ringing,
% etc should still be radially symmetric) with change in defocus, so this
% looks right.  At some point defocus will be large enough that you'll run
% into sampling problems.  You could check this by starting with very high
% sampling density and decreasing until just before the point where you
% start to see issues.

% Set up parameters structure
wvfParams0.zcoeffs = zeros(65,1);
wvfParams0.measpupilMM = 8;
wvfParams0.calcpupilMM = 3;
wvfParams0.wls = 550;
wvfParams0.nominalFocusWl = 550;
wvfParams0.defocusDiopters = 0;
wvfParams0.sizeOfFieldPixels = 201;
wvfParams0.sizeOfFieldMM = 16.212;
wavelengthOffset = 50;

nominalFocusWavelength = 550;
theWavelength = 550;
wavelengthOffset = 50;
sizeOfFieldPixels = 201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 3;
defocusDiopters = 0;

figure; clf; hold on
wvfParams1 = wvfParams0;
[wvfParams1] = wvfComputePSF(wvfParams0);
whichRow = floor(wvfParams1.sizeOfFieldPixels/2) + 1;
onedPSF1 = wvfParams1.psf(whichRow,:);
arcminutes = wvfParams1.arcminperpix*((1:wvfParams1.sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 2);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
wvfParams2 = wvfParams0;
wvfParams2.wls = wvfParams2.wls - wavelengthOffset;
[wvfParams2] = wvfComputePSF(wvfParams2);
onedPSF2 = wvfParams2.psf(whichRow,:);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
wvfParams3 = wvfParams0;
wvfParams3.wls = wvfParams3.wls + wavelengthOffset;
[wvfParams3] = wvfComputePSF(wvfParams3);
onedPSF3 = wvfParams3.psf(whichRow,:);
plot(arcminutes(index),onedPSF3(index),'g','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title(sprintf('Effect of defocus, +/- %0.1f nm',wavelengthOffset));

%% TEST3.  Include the Stiles-Crawford effect.
%
% Note from DHB.  I'm not sure here what the effect should look like, given
% that we're using a diffraction limited pupil, so this mostly just
% demonstrates that the code runs.
%
% Note from HH: This looks about right.
wvfParams0.zcoeffs = zeros(65,1);
wvfParams0.measpupilMM = 8;
wvfParams0.calcpupilMM = 8;
wvfParams0.wls = 550;
wvfParams0.nominalFocusWl = 550;
wvfParams0.defocusDiopters = 0;
wvfParams0.sizeOfFieldPixels = 201;
wvfParams0.sizeOfFieldMM = 16.212;
wvfParams0.sceParams = sceCreate(theWavelength,'berendshot');

figure; clf; hold on
wvfParams1 = wvfParams0;
wvfParams1.sceParams = [];
[wvfParams1] = wvfComputePSF(wvfParams1);
whichRow = floor(wvfParams1.sizeOfFieldPixels/2) + 1;
onedPSF1 = wvfParams1.psf(whichRow,:);
arcminutes = wvfParams1.arcminperpix*((1:wvfParams1.sizeOfFieldPixels)-whichRow);
index = find(abs(arcminutes) < 4);
plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
wvfParams2 = wvfParams0;
[wvfParams2] = wvfComputePSF(wvfParams2);
onedPSF2 = wvfParams2.psf(whichRow,:);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
title(sprintf('SCE, returned strehl %0.3f, direct from non-SCE diffrac %0.3f, frac. abs. is %0.2f',wvfParams2.strehl,max(wvfParams2.psf(:))/max(wvfParams1.psf(:)),wvfParams2.sceFrac));


%% TEST4.  Compute the psfs for the sample data
%
% We have a set of measured Zernike coefficients for 9 subjects.  Let's make sure
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
%
% Note from HH: Surprised by the the large changes with the 3mm pupil
% between zero defocus and that required to optimize the strehl.   Ran
% it with calcpupil of 6mm (since already had lots of calculations done
% for these subjects at this size) and with or without the SCE the
% values of defocus required to optimize the monochromatic strehl
% matched my calculations- at least within the resolution of your
% routine, so I guess the 3mm result is just what happens.  Not exactly
% and independent test, but at least verifies that these routines
% produce the same result as my original function. 
%
% Note from HH: For real calculations, using a defocus increment smaller
% than 0.25 Diopters would be wise.
wvfParams0.zcoeffs = zeros(65,1);
wvfParams0.measpupilMM = 6;
wvfParams0.calcpupilMM = 3;
wvfParams0.wls = 550;
wvfParams0.nominalFocusWl = 550;
wvfParams0.defocusDiopters = 0;
wvfParams0.sizeOfFieldPixels = 201;
wvfParams0.sizeOfFieldMM = 16.212;
theZernikeCoeffs = load('sampleZernikeCoeffs.txt');
DOSCE = 0;
if (DOSCE)
    wvfParams0.sceParams = sceCreate(theWavelength,'berendshot');
else
    wvfParams0.sceParams = sceCreate(theWavelength,'none');
end

figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
position(4) = 1600;
set(gcf,'Position',position);

% In this loop, 9 is number of subjects.  Hard coded so that plot comes out
% 3 by 3.  This counts on their being at least 9 subjects in the sample data
% file, which there are at present.
for i = 1:9
    subplot(3,3,i); hold on
    wvfParams1 = wvfParams0;
    wvfParams1.zcoeffs = zeros(61,1);
    wvfParams1 = wvfComputePSF(wvfParams1);
    whichRow = floor(wvfParams1.sizeOfFieldPixels/2) + 1;
    arcminutes = wvfParams1.arcminperpix*((1:wvfParams1.sizeOfFieldPixels)-whichRow);
    diffracPSF1 = wvfParams1.psf;
    onedPSF1 = wvfParams1.psf(whichRow,:);
    index = find(abs(arcminutes) < 6);
    plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
    
    wvfParams2 = wvfParams0;
    wvfParams2.zcoeffs = theZernikeCoeffs(:,i);
    wvfParams2 = wvfComputePSF(wvfParams2);
    thePSF2 = psfCenter(wvfParams2.psf);
    onedPSF2 = thePSF2(whichRow,:);
    plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
    
    thePSF2 = psfCenter(thePSF2);
    onedPSF2 = thePSF2(whichRow,:);
    strehlDirect2 = max(thePSF2(:))/max(diffracPSF1(:));
    
    % Optimize strehl by varying defocus
    bestStrehl = 0;
    defocusDiopters = -2:0.25:2;
    for j = 1:length(defocusDiopters)
        wvfParams3 = wvfParams0;
        wvfParams3.zcoeffs = theZernikeCoeffs(:,i);
        wvfParams3.defocusDiopters = defocusDiopters(j);
        wvfParams3 = wvfComputePSF(wvfParams3);
        if (wvfParams3.strehl > bestStrehl)
            bestStrehl = wvfParams3.strehl;
            bestPSF3 = psfCenter(wvfParams3.psf);
            bestDefocusDiopters = wvfParams3.defocusDiopters;
        end
    end
    onedPSF3 = bestPSF3(whichRow,:);
    plot(arcminutes(index),onedPSF3(index),'g','LineWidth',2);
    strehlDirect3 = max(bestPSF3(:))/max(diffracPSF1(:));
    
    % Optimize using function.  This optimizes mass within criterion
    % radius, not strehl, and uses parameter search not the
    % exhaustive search just above.  One could re-write the above
    % look to optimize the same thing as the search program, but
    % right now life just seems too short.
    wvfParams4 = wvfParams0;
    wvfParams4.zcoeffs = theZernikeCoeffs(:,i);
    wvfParams4.criterionFraction = 0.9;
    wvfParams4.optimizeWl = wvfParams4.wls(1);
    wvfParams4 = wvfComputeOptimizedPSF(wvfParams4);
    thePSF4 = psfCenter(wvfParams4.psf);
    onedPSF4 = thePSF4(whichRow,:);
    plot(arcminutes(index),onedPSF4(index),'k:','LineWidth',1);

    xlabel('Arc Minutes');
    ylabel('PSF');
    title(sprintf('Subject %d, strehl %0.2f (no defocus), %0.2f (defocus of %0.2f/%0.2f D)',i,wvfParams2.strehl,bestStrehl,bestDefocusDiopters,wvfParams4.defocusDiopters));
    drawnow;
    
    if (DOSCE)
        fprintf('Subject %i, with SCE correction\n',i);
    else
        fprintf('Subject %i, no SCE correction\n',i);
    end
    fprintf('\tNo defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect2,wvfParams2.strehl);
    fprintf('\tWith defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect3,bestStrehl);
    
    % Store results so we can play with them later.
    wvfParamsArray1(i) = wvfParams1;
    wvfParamsArray2(i) = wvfParams2;
    wvfParamsArray3(i) = wvfParams3;
        
end

%% TEST5.  Verify that PSF averaging function works correctly.  

% This looks about right by eye.
figure; clf; hold on
c = get(gcf,'Colormap');
psfsToAverage(:,:,1) = wvfParamsArray3(1).psf;
psfsToAverage(:,:,2) = wvfParamsArray3(8).psf;
averagePSF = psfAverageMultiple(psfsToAverage);
subplot(1,3,1);
mesh(psfsToAverage(:,:,1));
view([-16 10]);
xlim([0 201]); ylim([0 201]); zlim([0 3e-4]);
subplot(1,3,2);
mesh(psfsToAverage(:,:,2));
view([-16 10]);
xlim([0 201]); ylim([0 201]); zlim([0 3e-4]);
subplot(1,3,3);
mesh(averagePSF);
view([-16 10]);
xlim([0 201]); ylim([0 201]); zlim([0 3e-4]);






