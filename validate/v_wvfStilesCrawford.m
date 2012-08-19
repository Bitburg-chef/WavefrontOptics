%% v_wvfStilesCrawford
%
% Note from DHB.  I'm not sure here what the effect should look like, given
% that we're using a diffraction limited pupil, so this mostly just
% demonstrates that the code runs.
%
% Note from HH: This looks about right.

%%
s_initISET

% For plotting limits
maxMIN = 2;
maxMM  = 1;
waveIdx = 1;
theWavelength = 550;

%%
wvfParams = wvfCreate;
sceP = sceCreate(theWavelength,'berendschot_data','centered');

%% No Stiles Crawford effect
wvfParams = wvfSet(wvfParams,'sce params',[]);
wvfParams = wvfComputePSF(wvfParams);

vcNewGraphWin; 
hold on
wvfPlot(wvfParams,'1d psf angle','min',[],maxMIN);

% Make a graph of the PSF within 1 mm of center
% vcNewGraphWin;
% wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

%%  Include the SCE in place
wvfParams = wvfSet(wvfParams,'sce params',sceP);
wvfParams = wvfComputePSF(wvfParams);
% vcNewGraphWin;
% wvfPlot(wvfParams,'2dpsf space','mm',maxMM);

[f,p] = wvfPlot(wvfParams,'1d psf angle','min',[],maxMIN);
set(p,'color','b')
hold on

wList = wvfGet(wvfParams,'wave');
strehl = wvfGet(wvfParams,'strehl',wList);
title(sprintf('Strehl ratio:  %.1f',strehl));

% Figure this out
% title(sprintf('SCE, returned strehl %0.3f, direct from non-SCE diffrac %0.3f, frac. abs. is %0.2f',...
%     strehl,max(wvfParams2.psf(:))/max(wvfParams1.psf(:)),wvfParams2.sceFrac));

%% End
