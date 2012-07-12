%%  s_wvf2OI
%  Create an ISET optical image based on the point spread functions
%  calculated with the wavefront toolbox.
%
%  We are currently trying to figure out whether we have the units right.
%  It is surprising to see that the aberration from the wvf is so extreme.
%  It is much more extreme than the Marimont and Wandell 3mm pupil
%  analyses.  If!! we have it right here.
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize
s_initISET
maxUM = 10;

%% First build a diffraction limited case.
% Recall that we compare the diffraction limited case in
% v_wvfDiffractionPSF.m and found agreement between WVF, ISET, and PTB.
% Here we simply test whether we can convert a single wavelength
% pointspread in WVF into the same pointspread in ISET using the wvf2oi
% function.
wvfP = wvfCreate;
wvfP = wvfComputePSF(wvfP);
thisWave = wvfGet(wvfP,'wave');

% The conversion between minutes and 'um' is not right for the typical
% human focal length.  Look into this.
% wvfPlot(wvfP,'2d psf space','min',thisWave);
[u,p,f] = wvfPlot(wvfP,'2d psf space','um',thisWave);
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);

oiD = wvf2oi(wvfP,'shift invariant');
oiD = oiSet(oiD,'name','Diffraction limited');
vcAddAndSelectObject(oiD); oiWindow;
plotOI(oiD,'psf','um',thisWave)
set(gca,'xlim',[-maxUM maxUM],'ylim',[-maxUM maxUM]);


%% Now, try it with a few different wavelengths.
% Something wrong with multiple wavelength ...
% wave = (400:100:700); wave = wave(:);
% wvfP = wvfCreate('wave',wave,'name',sprintf('Diffraction-limited'));

%% Create a wavefront objects from the sample mean wvf data
% The data were collected by Thibos and are described in the wvfLoadHuman
% function and reference therein.
wave = (400:100:700); wave = wave(:);
pupilMM = 3; 

% First, try it for a diffraction limited system.
zCoefs = zeros(65,1);
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('DL-%d',pupilMM));
wvfP = wvfComputePSF(wvfP);
wvfPlot(wvfP,'2d psf space','min',500);
wvfPlot(wvfP,'2d psf space','um',500);

oiD  = wvf2oi(wvfP,'shift invariant');
oiD  = oiSet(oiD,'name','DL 3mm');

%
%zCoefs = wvfLoadHuman(pupilMM);
% wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
% oiD = wvf2oi(wvfP,'human');
% oiD = oiSet(oiD,'name','Human 3mm');

vcAddAndSelectObject(oiD); oiWindow;

%% Convert the wavefront structure to an ISET optical image (OI)
% You could save the OI if you like.  But for this purpose, we just leave
% it here.
oi = wvf2oi(wvfP,'human');
oi = oiSet(oi,'name','Human 3mm wvf');
%fname = fullfile(isetRootPath,'data','optics','wvfHuman30.mat');
%vcExportObject(oi,fname);

%% Create an ISET scene and run it through the oi we just created
% The wavefront measurements make a greenish image with almost no short
% wavelength contrast.
sweepS = sceneCreate('sweep');
sweepS = sceneSet(sweepS,'h fov',2);
oi = oiCompute(sweepS,oi);
vcAddAndSelectObject(oi); oiWindow;

%% The graph of the wavefront line spread at each wavelength
% This has almost no short wavelength contrast.  Seems impossible.
plotOI(oi,'ls wavelength');
title('Wavefront toolbox')

%% Repeat the process using the Marimont/Wandell estimate
% This is the one I included in my textbook.  It seems much more sensible
% to me.  So, we might be computing the wavefront wrong, possibly units.
% We are deep into it now.
oiMW = oiCreate('human');
oiMW = oiCompute(sweepS,oiMW);
vcAddAndSelectObject(oiMW); oiWindow;
plotOI(oiMW,'ls wavelength');
title('Marimont and Wandell')
