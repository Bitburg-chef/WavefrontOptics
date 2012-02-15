%% v_wvfSampleData 
% Compute psfs for the sample data
%
% We have Zernike coefficients measured in 9 subjects.  We compute the psfs
% using these, and look at a slice through each of them. 
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
% Note from HH: Surprised by the large changes with the 3mm pupil between
% zero defocus and that required to optimize the strehl.   Ran it with
% calcpupil of 6mm (since already had lots of calculations done for these
% subjects at this size) and with or without the SCE the values of defocus
% required to optimize the monochromatic strehl matched my calculations- at
% least within the resolution of your routine, so I guess the 3mm result is
% just what happens.  Not exactly an independent test, but at least
% verifies that these routines produce the same result as my original
% function. 
%
% Note from HH: For real calculations, using a defocus increment smaller
% than 0.25 Diopters would be wise.
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize

wvfParams0 = wvfCreate('measured pupil',6,'calculated pupil',3);

% Sample data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');

theZernikeCoeffs = load(sDataFile);

theZernikeCoeffs = theZernikeCoeffs(:,1:4);
nSubjects = size(theZernikeCoeffs,2);
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Stiles Crawford
DOSCE = 0;
theWavelength = 550;
if (DOSCE), wvfParams0.sceParams = sceCreate(theWavelength,'berendshot');
else        wvfParams0.sceParams = sceCreate(theWavelength,'none');
end

%%
vcNewGraphWin;
for ii = 1:nSubjects
    fprintf('** Subjects %d\n',ii)
    
    subplot(nRows,nCols,ii)
    
    % 
    wvfParams1 = wvfParams0;
    wvfParams1.zcoeffs = zeros(61,1);
    wvfParams1 = wvfComputePSF(wvfParams1);
    whichRow = floor(wvfParams1.sizeOfFieldPixels/2) + 1;
    arcminutes = wvfParams1.arcminperpix*((1:wvfParams1.sizeOfFieldPixels)-whichRow);
    diffracPSF1 = wvfParams1.psf;
    onedPSF1 = wvfParams1.psf(whichRow,:);
    index = find(abs(arcminutes) < 6);
    plot(arcminutes(index),onedPSF1(index),'r','LineWidth',4);
    hold on
    
    wvfParams2 = wvfParams0;
    wvfParams2.zcoeffs = theZernikeCoeffs(:,ii);
    wvfParams2 = wvfComputePSF(wvfParams2);
    thePSF2 = psfCenter(wvfParams2.psf);
    onedPSF2 = thePSF2(whichRow,:);
    plot(arcminutes(index),onedPSF2(index),'b','LineWidth',4);
    hold on
    
    thePSF2 = psfCenter(thePSF2);
    onedPSF2 = thePSF2(whichRow,:);
    strehlDirect2 = max(thePSF2(:))/max(diffracPSF1(:));
    
    % Optimize strehl by varying defocus
    bestStrehl = 0;
    defocusDiopters = -2:0.25:2;
    for j = 1:length(defocusDiopters)
        wvfParams3 = wvfParams0;
        wvfParams3.zcoeffs = theZernikeCoeffs(:,ii);
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
    wvfParams4.zcoeffs = theZernikeCoeffs(:,ii);
    wvfParams4.criterionFraction = 0.9;
    wvfParams4.optimizeWl = wvfParams4.wls(1);
    wvfParams4 = wvfComputeOptimizedPSF(wvfParams4);
    thePSF4 = psfCenter(wvfParams4.psf);
    onedPSF4 = thePSF4(whichRow,:);
    plot(arcminutes(index),onedPSF4(index),'k:','LineWidth',1);

    xlabel('Arc Minutes');
    ylabel('PSF');
    title(sprintf('Subject %d, strehl %0.2f (no defocus), %0.2f (defocus of %0.2f/%0.2f D)',ii,wvfParams2.strehl,bestStrehl,bestDefocusDiopters,wvfParams4.defocusDiopters));
    drawnow;
    
    if (DOSCE)
        fprintf('Subject %i, with SCE correction\n',ii);
    else
        fprintf('Subject %i, no SCE correction\n',ii);
    end
    fprintf('\tNo defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect2,wvfParams2.strehl);
    fprintf('\tWith defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect3,bestStrehl);
    
    % Store results so we can play with them later.
    wvfParamsArray1(ii) = wvfParams1;
    wvfParamsArray2(ii) = wvfParams2;
    wvfParamsArray3(ii) = wvfParams3;
        
end
