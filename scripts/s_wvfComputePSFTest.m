% s_wvfComputePSFTest
%
% OBSOLETE:  SEE v_wvf scripts
%
% We willl move the 2nd and 3rd part to validation scripts in the validate
% directory, too.
%
% Original notes:
%
% Tests the monochromatic PSFs using Zernike coefficients by comparing with
% PTB.  
%
% See also: wvfComputePSF, wvfComputePupilFunction,
%           sceCreate, wvfGetDefocusFromWavelengthDifference
%
% 3/15/12  MDL  Updated to use wvfSet and fieldSampleSizeMMperPixel
% 8/21/11  dhb  Wrote it, based on code provided by Heidi Hofer.
% 9/7/11   dhb  Got this working with wvfParams i/o.
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

%% Clear
s_initISET

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
onedPSF2 = AiryPattern(radians,wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1)));

%% Repeat the calculation with a wavelength offset

newWave = 400;  
wvfParams1 = wvfParams0;
wvfParams1 = wvfSet(wvfParams1,'wave',newWave);
wvfParams1 = wvfSet(wvfParams1,'in focus wavelength', newWave);
wvfParams = wvfComputePSF(wvfParams1);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN)
 
hold on
onedPSF2 = AiryPattern(radians,wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1)));

%% Repeat the calculation with a different pupil size
pupilMM = 7;   % In millimeters?
wvfParams2 = wvfParams0;
wvfParams2 = wvfSet(wvfParams2,'calculated pupil',pupilMM);
wvfParams = wvfComputePSF(wvfParams2);

vcNewGraphWin;
wvfPlot(wvfParams,'1d psf angle','min',maxMIN);
hold on;

onedPSF2 = AiryPattern(radians,wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1));
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',wvfGet(wvfParams,'calculated pupil'),wvfParams.wls(1)));

%% Put ISET diffraction comparisons here or insert above

%% END
%% Move this to a different script.  Move TEST3 to a different script, too. 
% TEST2: Change the nominal focus away from a specified wavelength
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

nominalFocusWavelength = 550;
theWavelength = 550;
fieldSampleSizeMMperPixel = 16.212/201;
sizeOfFieldMM = 16.212;
measpupilMM = 8;
calcpupilMM = 3;
defocusDiopters = 0;

% Set up parameters structure
wvfParams0 = wvfSet(wvfParams0,'zcoeffs',zeros(65,1));
wvfParams0 = wvfSet(wvfParams0,'measured pupil',measpupilMM);
wvfParams0 = wvfSet(wvfParams0,'calculated pupil',calcpupilMM);
wvfParams0 = wvfSet(wvfParams0,'wave',theWavelength);
wvfParams0 = wvfSet(wvfParams0,'infocus wavelength',nominalFocusWavelength);
wvfParams0 = wvfSet(wvfParams0,'defocus diopters',defocusDiopters);

% This looks like it might be redundant and should be removed
wvfParams0 = wvfSet(wvfParams0,'fieldsamplesize',fieldSampleSizeMMperPixel);
wvfParams0 = wvfSet(wvfParams0,'field size mm',sizeOfFieldMM);

wavelengthOffset = 50;
arcminPerPix = wvfGet(wvfParams1,'arcmin per pix');

figure; clf; hold on
wvfParams1 = wvfComputePSF(wvfParams0);
whichRow = floor(wvfGet(wvfParams1,'npixels')/2) + 1;
onedPSF1 = wvfParams1.psf(whichRow,:);
arcminutes = arcminPerPix*((1:wvfGet(wvfParams1,'npixels'))-whichRow);
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
wvfParams0.fieldSampleSizeMMperPixel = 16.212/201;
wvfParams0.sizeOfFieldMM = 16.212;
wvfParams0.sceParams = sceCreate(theWavelength,'berendschot_data');

figure; clf; hold on
wvfParams1 = wvfParams0;
wvfParams1 = rmfield(wvfParams1,'sceParams');
% wvfParams1.sceParams = [];
[wvfParams1] = wvfComputePSF(wvfParams1);
whichRow = floor(wvfGet(wvfParams1,'npixels')/2) + 1;
onedPSF1 = wvfParams1.psf(whichRow,:);
arcminutes = arcminPerPix*((1:wvfGet(wvfParams1,'npixels'))-whichRow);
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
wvfParams0.fieldSampleSizeMMperPixel = 16.212/201;
wvfParams0.sizeOfFieldMM = 16.212;
theZernikeCoeffs = importdata('sampleZernikeCoeffs.txt');
DOSCE = 0;
if (DOSCE)
    wvfParams0.sceParams = sceCreate(theWavelength,'berendschot_data');
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
    whichRow = floor(wvfGet(wvfParams1,'npixels')/2) + 1;
    arcminutes = arcminPerPix*((1:wvfGet(wvfParams1,'npixels'))-whichRow);
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






