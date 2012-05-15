%% v_wvfSampleData
%
% Compute psfs for the sample data
%
% HH provided Zernike coefficients measured in 9 subjects.  We compute the
% psfs using these, and look at a slice through each of them.  We try
% things various different ways.
%
% The computed PSFs are recentered, with their maximum in the center, so that we
% see the real peak when we take the 1D slice.
%
% This also optimizes the defocus to maximize the strehl ratio for each
% subject, so you can see the (large) effect of doing that.
%
% Note from DHB.  Again, I don't know if these are correct, but at least
% you can see that you get a wide range of answers by using different
% subjects' data.
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
s_initISET
waveIdx = 1;
maxMIN = 6;

%% Use Heidi Hofer's sample data here

% Set values in millimeters
wvfP = wvfCreate('measured pupil',6,'calculated pupil',3);

% Sample data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = load(sDataFile);

whichSubjects = 1:2; nSubjects = length(whichSubjects);
theZernikeCoeffs = theZernikeCoeffs(:,whichSubjects);

% For plotting
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Stiles Crawford
DOSCE = 0;
sceWavelength = 550;
if (DOSCE), wvfP.sceParams = sceCreate(sceWavelength,'berendshot');
else        wvfP.sceParams = sceCreate(sceWavelength,'none');
end

%%
vcNewGraphWin;
for ii = 1:nSubjects
    fprintf('** Subject %d\n',ii)
    
    % Compute the diffraction limited version of the PSF
    wvfP = wvfSet(wvfP,'zcoeffs',zeros(61,1));
    wvfP = wvfComputePSF(wvfP);
    % Diffraction limited
    udataD = wvfPlot(wvfP,'1d psf angle','min',waveIdx,maxMIN);
    hold on;
    
    % Now, set it up for the typical subject
    wvfP = wvfSet(wvfP,'zcoeffs',theZernikeCoeffs(:,ii));
    wvfP = wvfComputePSF(wvfP);
    
    [udataS, pData] = wvfPlot(wvfP,'1d psf angle','min',waveIdx,maxMIN);
    set(pData,'color','b');
    hold on;
    
    strehlDirect = max(udataS.y(:))/max(udataD.y(:));
    fprintf('Strehl ratio with no defocus:  %.3f\n',strehlDirect);
    
    
end

%% End


% Optimize strehl by varying defocus
%     bestStrehl = 0;
%     defocusDiopters = -2:0.25:2;
%     thisStrehl = zeros(1,length(defocusDiopters));
%     for jj = 1:length(defocusDiopters)
%         %         wvfParams3 = wvfP;
%         %         wvfParams3.zcoeffs = theZernikeCoeffs(:,ii);
%         wvfParams = wvfSet(wvfParams,'defocus diopters',defocusDiopters(jj));
%         wvfParams = wvfComputePSF(wvfParams);
%         thisStrehl(jj) = wvfGet(wvfParams,'strehl');
%         if thisStrehl(jj) == max(thisStrehl)
%             wvfParamsBest = wvfParams;
%         end
%     end

% Show the best optical correction (defocus) for this pointspread
%     vcNewGraphWin; plot(defocusDiopters,thisStrehl,'-o');
%     title('Defocus effect on peak of the PSF');

% Best defocus
%     figure(f); subplot(nRows,nCols,ii)
%     wvfPlot(wvfParamsBest,'2d psf angle','min',maxMIN);
%
% A little slow, probably interesting.  Deal with this later.
% Perhaps put it in a separate script.
% - BW
% Optimize using function.  This optimizes mass within criterion
% radius, not strehl, and uses parameter search not the
% exhaustive search just above.  One could re-write the above
% look to optimize the same thing as the search program, but
%     % right now life just seems too short.
%     wvfParams4 = wvfP;
%     wvfParams4.zcoeffs = theZernikeCoeffs(:,ii);
%     wvfParams4.criterionFraction = 0.9;
%     wvfParams4.optimizeWl = wvfParams4.wls(1);
%     wvfParams4 = wvfComputeOptimizedPSF(wvfParams4);
%     thePSF4 = psfCenter(wvfParams4.psf);
%     onedPSF4 = thePSF4(whichRow,:);
%     plot(arcminutes(index),onedPSF4(index),'k:','LineWidth',1);
%
%     xlabel('Arc Minutes');
%     ylabel('PSF');
%     title(sprintf('Subject %d, strehl %0.2f (no defocus), %0.2f (defocus of %0.2f/%0.2f D)',ii,wvfParams2.strehl,bestStrehl,bestDefocusDiopters,wvfParams4.defocusDiopters));
%     drawnow;
%
%     if (DOSCE)
%         fprintf('Subject %i, with SCE correction\n',ii);
%     else
%         fprintf('Subject %i, no SCE correction\n',ii);
%     end
%     fprintf('\tNo defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect2,wvfParams2.strehl);
%     fprintf('\tWith defocus: direct strehl %0.3f, returned, %0.3f\n',strehlDirect3,bestStrehl);
%
%     % Store results so we can play with them later.
%     wvfParamsArray1(ii) = wvfParams1;
%     wvfParamsArray2(ii) = wvfParams2;
%     wvfParamsArray3(ii) = wvfParams3;
