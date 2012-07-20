% v_wvfDiffractionPSF
%
% Tests the monochromatic PSFs computed from Zernike coefficients. The
% script compares the computations with those in PTB.
%
% See also: wvfCreate, wvfGet, wvfSet, wvfComputePSF, wvfComputePupilFunction,
%   wvfLCAFromWavelengthDifference
%
% 6/5/12  dhb  Got this to work with current code.
% 7/1/12  bw   Replaced waveIdx with wList in wvfGet and here and ...
%              Maybe we should just use GIT for these comments.
% 7/4/12  dhb  Remove some unnecessary lines (don't need to set PSF_STALE,
%              for example.)
% 7/16/12 dhb  Remove dependency on rad2deg to avoid namespace conflicts.
%         dhb  Get rid of extraneous calls to make new plotting windows.
%              These were producing blank figures when the script ran,
%              because the wvfPlot function creates its own windows.
%         dhb  Verify that diffraction limited output is independent of
%              specified measured pupil size.
%         dhb  Only compute VSET version if it is on path.  This is
%              a bit of a tilt at windmills, since there are many other
%              vset dependencies.
% 7/20/12 dhb  Change iset comments to vset.
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize
s = which('v_wvfDiffractionPSF');
cd(fileparts(s));
clear; close all;
s_initISET;

%% Compare pointspread function in wvf with psf in Psych Toolbox

% When the Zernike coefficients are all zero, the wvfComputePSF code should
% return the diffraction limited PSF.  We test whether this works by
% comparing to the diffraction limited PSF implemented in the PTB routine
% AiryPattern.

% Set up parameters structure
measPupilMM = 4;
calcPupilMM = 3;
wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'measured pupil size',measPupilMM);
wvf0 = wvfSet(wvf0,'calc pupil size',calcPupilMM);

% Plotting ranges for MM, UM, and Minutes of angle
maxMM = 1;
maxUM = 20;
maxMIN = 2;

% Which wavelength index (wave(idx)) (there is only one) to plot
wList = wvfGet(wvf0,'wave');

% Calculate the PSF, normalized to peak of 1.
wvf0 = wvfComputePSF(wvf0);

% Make sure psf computed this way (with zcoeffs zeroed) matches
% what is returned by our internal get of diffraction limited
% psf.
psf = wvfGet(wvf0,'psf');
diffpsf = wvfGet(wvf0,'diffraction psf');
if (any(abs(psf(:)-diffpsf(:)) ~= 0))
    error('Internal computation of diffraction limited psf does not match explicit calc');
end

% Get out parameters for various checks
calcWavelength = wvfGet(wvf0,'wavelength');
measWavelength = wvfGet(wvf0,'measured wavelength');
if (measWavelength ~= calcWavelength)
    error('Measured and calculation wavelengths should match at this point');
end

% Make a graph of the PSF within maxUM of center
wvfPlot(wvf0,'2dpsf space','um',wList,maxUM);

% Make a graph of the PSF within 2 arc min
wvfPlot(wvf0,'2dpsf angle','min',wList,maxMIN);

% Plot the middle row of the psf, scaled to peak of 1
wvfPlot(wvf0,'1d psf angle normalized','min',wList,maxMIN);
hold on

% Used for plotting comparisons below
arcminutes = wvfGet(wvf0,'support arcmin','min',wList);
arcminpersample = wvfGet(wvf0,'ref psf sample interval');
arcminpersample1 = wvfGet(wvf0,'psf arcmin per sample',wList);
arcminpersample2 = wvfGet(wvf0,'psf angle per sample',[],wList);
if (arcminpersample1 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end
if (arcminpersample2 ~= arcminpersample1)
    error('Default units of get on ''psfanglepersample'' unexpectedly changed');
end
index = find(abs(arcminutes) < 2);
radians = (pi/180)*(arcminutes/60);

% Compare to what we get from PTB AiryPattern function -- should match
onedPSF2 = AiryPattern(radians,calcPupilMM ,calcWavelength);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
figNum = gcf;
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',calcPupilMM,calcWavelength));

%% Do the same thing using vset functions, if they exist on the path
%
% The conversion between VSET and these other methods is pretty good, too.
% So, diffraction limited point spread, measured for 3.0mm is the same when
% done with VSET, PTB and WVF.
if (exist('oiCreate','file'))
    thisWave = 550;
    oi = oiCreate;
    optics = oiGet(oi,'optics');
    fLength = 0.017;              % Human flength is about 17 mm
    fNumber = 17/calcPupilMM;     % Set pupil diameter
    
    optics = opticsSet(optics,'flength',fLength);  % Roughly human
    optics = opticsSet(optics,'fnumber',fNumber);   % Roughly human
    oi = oiSet(oi,'optics',optics);
    [uData,g] = plotOI(oi,'psf',[],thisWave); close(g);
    
    figure(figNum)
    [r,c] = size(uData.x);
    mid = ceil(r/2);
    psfMid = uData.psf(mid,:);
    posMM = uData.x(mid,:)/1000;               % Microns to mm
    posMinutes = 60*(180/pi)*(atan2(posMM,opticsGet(optics,'flength','mm')));
    
    g = wvfPlot(wvf0,'1d psf angle normalized','min',wList,maxMIN);
    hold on
    plot(posMinutes,psfMid/max(psfMid(:)),'ko')
    hold on
    plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
    xlabel('Arc min')
    set(gca,'xlim',[-2 2])
    grid on
    legend('WVF','VSET','PTB');
end

%% Repeat the calculation with a wavelength offset.  
% To keep the new wavelength in focus in the calculations, we add an
% explicit observer focus correction, with the amount computed by
% wvfLCAFromWavelengthDifference relative to the measured wavelength
wList = 400;  
wvf1 = wvf0;
wvf1 = wvfSet(wvf1,'wave',wList);
lcaDiopters = wvfLCAFromWavelengthDifference(wList,wvfGet(wvf1,'measured wl'));
wvf1 = wvfSet(wvf1,'calc observer focus correction',lcaDiopters);
wvf1 = wvfComputePSF(wvf1);

w = wvfGet(wvf1,'wave');
pupilSize = wvfGet(wvf1,'calcpupilsize','mm');

wvfPlot(wvf1,'1d psf angle normalized','min',w,maxMIN);
hold on
onedPSF2 = AiryPattern(radians,pupilSize,w);
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,w));

% PSF angular sampling should be the same across wavelengths
arcminpersample2 = wvfGet(wvf1,'psf angle per sample','min',w);
if (arcminpersample2 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end

%% Repeat the calculation with a different pupil size at original wavelength
pupilMM = 7; 
wvf2  = wvf0;
wvf2 = wvfSet(wvf2,'measured pupil size',pupilMM);
wvf2  = wvfSet(wvf2,'calc pupil size',pupilMM);
wvf2  = wvfComputePSF(wvf2);
wList = wvfGet(wvf2,'wave');
pupilSize = wvfGet(wvf2,'calc pupil size','mm');

% Compare the curves
wvfPlot(wvf2,'1d psf angle normalized','min',wList,maxMIN);
onedPSF2 = AiryPattern(radians,pupilSize,wList);

hold on
plot(arcminutes(index),onedPSF2(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,wList));



