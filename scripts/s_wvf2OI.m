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

%% Create a wavefront objects from the sample mean wvf data
% The data were collected by Thibos and are described in the wvfLoadHuman
% function and reference therein.
pupilMM = 3; zCoefs = wvfLoadHuman(pupilMM);
wave = (400:10:700); wave = wave(:);
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));

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
