%% s_wvfHumanScenes


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
