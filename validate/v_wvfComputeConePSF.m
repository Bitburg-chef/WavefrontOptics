% v_wvfComputeConePSF
%
% Test the routines that compute L, M, and S cone PSFs from Zernike
% coefficients.
%
% The diffraction limited calcs seem to match up with their
% Figure 2 pretty well, the MTFs agree well, and Figure 4B
% is reproduced to good approximation.  Still need to check
% Figure 4a.
%
% See also: wvfComputeConePSF, wvfComputePSF, wvfComputePupilFunction,
%   sceGetParams, wvfGetDefocusFromWavelengthDifference
%
% The circular averaging is not a good idea for a single subject, but if
% you want to obtain an average over subjects it might be good.
%
% 8/21/11  dhb  Wrote it.
% 3/15/12  mdl  Edited to use wvfSet. Also updated to use fieldSampleSizeMMperPixel
% 7/20/12  dhb  Got TEST1 to work without crashing, and possibly even to be correct.
% 7/23/12  dhb  OTF plot is looking vaguely reasonable.
%               Added Autrusseau equal energy OTFs for comparison

%% Initialize
clear; close all;
s = which('v_wvfComputeConePSF');
cd(fileparts(s));
s_initISET;

%% Parameters
%
% Did Autruesseau et al. incorporate a model of the SCE?
DOSCE = 0;
CIRCULARLYAVERAGE = 0;
CENTER = 0;
plotLimit = 6;
plotLimitFreq = 80;

%% Test data
dataSource = 'AutrusseauStandard';
switch (dataSource)
    case 'AutrusseauStandard';
        % This is the Autrussea standard observer.
        % I think their first coefficient is the
        % j = 0 Zernike mode number, and we don't use
        % this.  So we lop it off.
        %
        % Their coefficients are for a measured
        % pupil of 6 mm and that their calculations are
        % also for 6 mm. They use 570 nm as their
        % measured (in focus) wavelength.  [Their methods,
        % p. 2284].
        whichSubject = 1;
        dataFile = 'autrusseauStandardObserver.txt';
        theZernikeCoeffs = importdata(dataFile);
        theZernikeCoeffs = theZernikeCoeffs(2:15);
        %theZernikeCoeffs(4) = 0;
        measPupilMM = 6;
        calcPupilMM = 6;
        measWavelength = 570;
    case 'ThibosStatisticalModelMean'
        % This is the mean data from the Thibos model
        % for a 6 mm pupil.  The methods of the Autrusseau
        % indicate that they only used the first 15 (starting
        % at j = 0 coefficients.
        whichSubject = 1;
        load('IASstats60','sample_mean');
        theZernikeCoeffs = sample_mean;
        theZernikeCoeffs = theZernikeCoeffs(2:15);
        %theZernikeCoeffs(4) = 0;
        measPupilMM = 6;
        calcPupilMM = 6;
        measWavelength = 570;
end

% Cone sensitivities and equal energy weighting spectrum
load('T_cones_ss2');
conePsfInfo.S = S_cones_ss2;
conePsfInfo.T = T_cones_ss2;
conePsfInfo.spdWeighting = ones(conePsfInfo.S(3),1);

% Calculation wavelengths for PSF.  
wls = SToWls([400 10 31]);

%% TEST1: Try to reproduce Autrusseau et al. results.
% No focus optimization.
wvf0 = wvfCreate;

% Set important parameters
wvf0 = wvfSet(wvf0,'measured pupil size',measPupilMM);
wvf0 = wvfSet(wvf0,'calc pupil size',calcPupilMM);
wvf0 = wvfSet(wvf0,'zcoeffs',theZernikeCoeffs(:,whichSubject));
wvf0 = wvfSet(wvf0,'measured wavelength',measWavelength);
wvf0 = wvfSet(wvf0,'calc wavelengths',wls);
wvf0 = wvfSet(wvf0,'calc cone psf info',conePsfInfo);

% The short wavelength PSFs get blurred enough that we need more pixels to
% contain them than our defaults provide.  Adjust to keep sempling density
% in psf plane the same, but increase space sampled.
origSamples = wvfGet(wvf0,'spatial samples');
newSamples = 601;
wvf0 = wvfSet(wvf0,'spatial samples',newSamples);
%wvf0 = wvfSet(wvf0,'ref pupil plane size',newSamples/origSamples*wvfGet(wvf0,'ref pupil plane size'));
fprintf('Sampling pupil plane/psf with %d pixels\n',wvfGet(wvf0,'spatial samples'));
fprintf('Pupil plane info\n');
for wavelength = [400 500 600 700];
    fprintf('\t%d nm, %0.1f mm, %0.3f mm/pixel\n',...
        wavelength,wvfGet(wvf0,'pupil plane size','mm',wavelength),wvfGet(wvf0,'pupil plane size','mm',wavelength)/wvfGet(wvf0,'spatial samples'));
end
fprintf('PSF plane info\n');
for wavelength = [400 500 600 700];
    fprintf('\t%d nm, %0.1f minutes, %0.3f min/pixel\n',...
        wavelength,wvfGet(wvf0,'psf angle per sample','min',wavelength)*wvfGet(wvf0,'spatial samples'),wvfGet(wvf0,'psf angle per sample','min',wavelength));
end

if (DOSCE == 1)
    sce = sceCreate(wls,'berendshot');
    wvf0 = wvfSet(wvf0,'sce params',sce);
else
    sce = sceCreate(wls,'none');
    wvf0 = wvfSet(wvf0,'sce params',sce);
end

% Compute LMS psfs both for a subject and diffraction limited
wvfParams1 = wvf0;
wvfParams1 = wvfComputePSF(wvfParams1);
conePsf1 = wvfGet(wvfParams1,'cone psf');

wvfParams2 = wvf0;
wvfParams2 = wvfSet(wvfParams2,'zcoeffs',zeros(65,1));
wvfParams2 = wvfComputePSF(wvfParams2);
conePsf2 = wvfGet(wvfParams2,'cone psf');

if (CENTER)
    lpsf = psfCenter(conePsf1(:,:,1));
    mpsf = psfCenter(conePsf1(:,:,2));
    spsf = psfCenter(conePsf1(:,:,3));
    lpsfd = psfCenter(conePsf2(:,:,1));
    mpsfd = psfCenter(conePsf2(:,:,2));
    spsfd = psfCenter(conePsf2(:,:,3));
else
    lpsf = conePsf1(:,:,1);
    mpsf = conePsf1(:,:,2);
    spsf = conePsf1(:,:,3);
    lpsfd = conePsf2(:,:,1);
    mpsfd = conePsf2(:,:,2);
    spsfd = conePsf2(:,:,3);
end

if (CIRCULARLYAVERAGE)
    lpsf = psfCircularlyAverage(lpsf);
    mpsf = psfCircularlyAverage(mpsf);
    spsf = psfCircularlyAverage(spsf);
    lpsfd = psfCircularlyAverage(lpsfd);
    mpsfd = psfCircularlyAverage(mpsfd);
    spsfd = psfCircularlyAverage(spsfd);
end

whichRow = wvfGet(wvfParams1,'middle row');
for i = 1:length(wls)
    if (wvfGet(wvfParams1,'psf arcmin per sample',wls(1)) ~= wvfGet(wvfParams1,'psf arcmin per sample',wls(i)))
        error('Error in spatial sampling consistency across wavelengths');
    end
end
arcminutes = wvfGet(wvfParams1,'psf arcmin per sample',wls(1))*((1:wvfGet(wvfParams1,'spatial samples'))-whichRow);

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

%% Take a look in frequency domain.
%
% I started using psf2otf, but it 
% centers its output in a way I find
% counterintuitive so I just went back 
% to the straight fft.
lotf = fftshift(fft2(lpsf));
motf = fftshift(fft2(mpsf));
sotf = fftshift(fft2(spsf));
lotfd = fftshift(fft2(lpsfd));
motfd = fftshift(fft2(mpsfd));
sotfd = fftshift(fft2(spsfd));

% Figure out the scale in the frequency domain.
% THIS NEEDS DOUBLE CHECKING.
totalDegrees = (arcminutes(end)-arcminutes(1))/60;
cyclesDegreePerPixel = 1/totalDegrees;
cyclesdegree = cyclesDegreePerPixel*((1:wvfGet(wvfParams1,'spatial samples'))-whichRow);

%% Read in Autrusseau data for comparison
%
% The fields come in badly labeled, because I misunderstood
% the graph when I digitzed it.  The mtf in the file is not
% log mtf, it's just the straight mtf.
%
% Note that Autrussea et al. didn't use the Fourier transform
% instead literally convolved the psf at each wavelength
% with a sinusoidal stimulus at the same wavelength, and then summed
% up the results over wavelength, weighting by the cone sensitivities.
% I think what we do is equivalent, but we may have to code it up and
% check.
autrusseauFigure11 = ReadStructsFromText('autrusseauFigure11.txt');

% Make a plot of the horizontal MTF
theFig = figure; clf;
position = get(gcf,'Position');
position(3) = 1600;
set(gcf,'Position',position);
subplot(1,3,1); hold on
onedLOTFH = abs(lotf(whichRow,:));
onedLOTFV = abs(lotf(:,whichRow));
onedLOTFD = abs(lotfd(whichRow,:));
index = find(abs(cyclesdegree) < plotLimitFreq);
plot(cyclesdegree(index),log10(onedLOTFH(index)),'r','LineWidth',2);
plot(cyclesdegree(index),log10(onedLOTFV(index)),'r:','LineWidth',2);
plot([autrusseauFigure11.sf_cpd],log10([autrusseauFigure11.log10_Lmtf_ees]),'ro','MarkerSize',8,'MarkerFaceColor','r');
plot(cyclesdegree(index),log10(onedLOTFD(index)),'k','LineWidth',1);
xlim([0 plotLimitFreq]);
ylim([-3 0]);
xlabel('Cycles/Degree');
ylabel('LOG10 OTF');
if (CIRCULARLYAVERAGE)
    title('Circularized L cone PSF');
else
    title('L cone OTF');
end

subplot(1,3,2); hold on
onedMOTFH = abs(motf(whichRow,:));
onedMOTFV = abs(motf(:,whichRow));
onedMOTFD = abs(motfd(whichRow,:));
index = find(abs(cyclesdegree) < plotLimitFreq);
plot(cyclesdegree(index),log10(onedMOTFH(index)),'g','LineWidth',2);
plot(cyclesdegree(index),log10(onedMOTFV(index)),'g:','LineWidth',2);
plot([autrusseauFigure11.sf_cpd],log10([autrusseauFigure11.log10_Mmtf_ees]),'go','MarkerSize',8,'MarkerFaceColor','g');
plot(cyclesdegree(index),log10(onedMOTFD(index)),'k','LineWidth',1);
xlim([0 plotLimitFreq]);
ylim([-3 0]);
xlabel('Cycles/Degree');
ylabel('LOG10 OTF');
if (CIRCULARLYAVERAGE)
    title('Circularized M cone PSF');
else
    title('M cone OTF');
end

subplot(1,3,3); hold on
onedSOTFH = abs(sotf(whichRow,:));
onedSOTFV = abs(sotf(:,whichRow));
onedSOTFD = abs(sotfd(whichRow,:));
index = find(abs(cyclesdegree) < plotLimitFreq);
plot(cyclesdegree(index),log10(onedSOTFH(index)),'b','LineWidth',2);
plot(cyclesdegree(index),log10(onedSOTFV(index)),'b:','LineWidth',2);
plot([autrusseauFigure11.sf_cpd],log10([autrusseauFigure11.log10_Smtf_ees]),'bo','MarkerSize',8,'MarkerFaceColor','b');
plot(cyclesdegree(index),log10(onedSOTFD(index)),'k','LineWidth',1);
xlim([0 plotLimitFreq]);
ylim([-3 0]);
xlabel('Cycles/Degree');
ylabel('LOG10 OTF');
if (CIRCULARLYAVERAGE)
    title('Circularized S cone PSF');
else
    title('S cone OTF');
end
drawnow;

%% Make a figure comparable to Autrusseau et al, Figure 2
% (top row) and Figure 4b (bottom row).
%
% This shows diffraction limited and standard observer PSFs at different
% wavelengths, given focus at 570.

wavelengths = [400 550 700];
figure; clf;
for i = 1:length(wavelengths);
    wavelength = wavelengths(i);
    
    subplot(2,length(wavelengths),i); hold on
    [nil,p] = wvfPlot(wvfParams2,'image pupil phase','mm',wavelength,'no window');
    
    focusWl = wvfGet(wvfParams2,'measured wavelength');
    subplot(2,length(wavelengths),i+length(wavelengths)); hold on
    psf = wvfGet(wvfParams2,'psf',wavelength);
    maxVal = max(psf(:));
    % [nil,p] = wvfPlot(wvfParams2,'2d psf angle','min',wavelength,'no window');
    [nil,p] = wvfPlot(wvfParams2,'image psf angle','min',wavelength,'no window');
    h = get(p,'Parent');
    view([0 90]); ylim([-30 30]); xlim([-30 30]); axis('square');
    title(sprintf('%d nm, focus %d nm, max = %0.5f',wavelength,focusWl,maxVal));
end

figure; clf;
for i = 1:length(wavelengths);
    wavelength = wavelengths(i);
    
    subplot(2,length(wavelengths),i); hold on
    [nil,p] = wvfPlot(wvfParams1,'image pupil phase','mm',wavelength,'no window');
    %h = get(p,'Parent');
    %view([0 90]); ylim([-30 30]); xlim([-30 30]); axis('square');
    %title(sprintf('%d nm, focus %d nm, max = %0.5f',wavelength,focusWl,maxVal));

    subplot(2,length(wavelengths),i+length(wavelengths)); hold on
    psf = wvfGet(wvfParams1,'psf',wavelength);
    maxVal = max(psf(:));
    %[nil,p] = wvfPlot(wvfParams1,'2d psf angle','min',wavelength,'no window');
    [nil,p] = wvfPlot(wvfParams1,'image psf angle','min',wavelength,'no window');
    h = get(p,'Parent');
    view([0 90]); ylim([-30 30]); xlim([-30 30]); axis('square');
    title(sprintf('%d nm, focus %d nm, max = %0.5f',wavelength,focusWl,maxVal));
end

return

%% TEST2.  Optimize focus and add to the plot.
%
% This takes a long time.


%Should be using sets/gets

wvfParams3 = wvf0;
wvfParams3.coneWeights = [1 1 0];
wvfParams3.criterionFraction = 0.9;

% This takes a long time and produces an error that could be fixed by BW,
% but he is too lazy.
%  Error using wvfGet (line 590)
%   Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF

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





