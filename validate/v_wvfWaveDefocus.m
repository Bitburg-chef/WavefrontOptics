% v_wvfWaveDefocus
%
% Compute for wavelengths different from nominal focus wavelength
%
% The psf should get broader away from the focus wavelength. How much
% broader, I don't know but we can at least verify the qualitative
% behavior, again for the diffaction limited case.
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
%
% See also:  v_wvfDiffractionPSF, v_wvfWaveDefocus
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize and set parameters

s_initISET

% Ranges for plotting
maxMIN = 2;
maxMM  = 1;
maxUM  = 20;
waveIdx = 1;

%% Calculate point spread for wavelength defocus

% Set up default parameters structure with diffraction limited default
wvfP = wvfCreate;
wave = 400:50:700;
wvfP = wvfSet(wvfP,'wave',wave);
nWave = wvfGet(wvfP,'n wave');

%% Compute and plot the default
wvfParams = wvfComputePSF(wvfP);

%% Plot the series of lines
f = vcNewGraphWin([],'tall');
for ii=1:nWave
    subplot(nWave,1,ii)
    [~,p] = wvfPlot(wvfParams,'1d psf space','um',ii,maxUM);
    title(sprintf('wave %d',wave(ii)));
end

% Alternative plotting method
% [~,p] = wvfPlot(wvfParams,'1d psf angle','min',ii,maxMIN);

%% End

% Used to do the stuff below.  Then did a major re-write.  Delete at some
% point. 
%
% Offset the stimulus wavelength from the in focus wavelength
% inFocus = wvfGet(wvfP,'inFocus Wavelength');
% wavelengthOffset = 50;
% wvfParams = wvfSet(wvfP,'wave',inFocus - wavelengthOffset);
% wvfParams = wvfComputePSF(wvfParams);
% 
% vcNewGraphWin(f); 
% hold on
% [~,p] = wvfPlot(wvfParams,'1d psf angle','min',waveIdx,maxMIN);
% set(p,'color','b','linestyle','--')
% str{2} = sprintf('wave %d',wvfGet(wvfParams,'wave'));
% 
% % vcNewGraphWin;
% % wvfPlot(wvfParams,'2dpsf space','mm',maxMM);
% 
% %% Now shift the offset the other way
% 
% inFocus = wvfGet(wvfP,'inFocus Wavelength');
% wavelengthOffset = 50;
% wvfParams = wvfSet(wvfP,'wave',inFocus + wavelengthOffset);
% wvfParams = wvfComputePSF(wvfParams);
% 
% hold on
% [~,p] = wvfPlot(wvfParams,'1d psf angle','min',waveIdx,maxMIN);
% set(p,'color','r','linestyle',':')
% str{3} = sprintf('wave %d',wvfGet(wvfParams,'wave'));
% legend(str)
% title(sprintf('Effect of defocus, + %0.1f nm',wavelengthOffset));

% vcNewGraphWin;
% wvfPlot(wvfParams,'2dpsf space','mm',maxMM);


