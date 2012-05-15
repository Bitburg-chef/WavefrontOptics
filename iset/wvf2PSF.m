function [psf, wave, umPerSample, wvfP] = wvf2PSF(wvfP)
% Convert a wvf structure to ISET shift-invariant PSF data
%
%    [psf, wave, umPerSample, wvfP] = wvf2PSF(wvfP)
%
% For each wavelength in wvf, compute the PSF with proper units and
% place it in an ISET shift-invariant PSF format that can be used for human
% optics simulation.
%
% The wvfP is a standard wavefront optics toolbox structure.
% The psf is computed at the wave values in the structure.  The updated
% structure with the PSFs can be returned.
%
% The ISET data can be saved using ieSaveSIDataFile as in the example
% below, which loads the standard human data for a particular pupil size.
%
% Example:
%    pupilMM = 3; zCoefs = wvfLoadHuman(pupilMM);
%    wave = [450:100:650]';
%    wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
%
%    [psf, wave, umPerSample, wvfP] = wvf2PSF(wvfP);
%    fName = sprintf('psfSI-%s',wvfGet(wvfP,'name'));
%    ieSaveSIDataFile(psf,wave,umPerSample,fName);
%
%    vcNewGraphWin([],'tall');
%    subplot(2,1,1), wvfPlot(wvfP,'image psf','um',1,30);
%    subplot(2,1,2), wvfPlot(wvfP,'image psf','um',2,30);
%
% Copyright Imageval 2012

%% Parameters
if ieNotDefined('wvfP'), error('wvf parameters required.'); end
wave = wvfGet(wvfP,'wave');
nWave = wvfGet(wvfP,'nwave');

% Use WVF to compute the PSFs
wvfP = wvfComputePSF(wvfP);

% Set up to interpolate the PSFs for ISET
% number of pixels and spacing in microns between samples
nPix = 256;                 % Number of pixels
umPerSample = [0.25,0.25];  % In microns
iSamp = (1:nPix)*umPerSample(1);
iSamp = iSamp - mean(iSamp);
iSamp = iSamp(:);
psf = zeros(nPix,nPix,nWave);

for ii=1:nWave,
    thisPSF = wvfGet(wvfP,'psf',ii);  % vcNewGraphWin; imagesc(thisPSF)
    samp = wvfGet(wvfP,'samples space','um',ii);
    samp = samp(:);
    psf(:,:,ii) = interp2(samp,samp',thisPSF,iSamp,iSamp');
    % wvfPlot(wvfP,'image psf space','um',ii,50)
end     


end


%% From ISET script on how to create an SI data file

% Now, write out a file containing the relevant point spread function
% data, along with related variables.
% umPerSample = [0.25,0.25];                % Sample spacing
% 
% % Point spread is a little square in the middle of the image
% h = zeros(128,128); h(48:79,48:79) = 1; h = h/sum(h(:));
% for ii=1:length(wave), psf(:,:,ii) = h; end     % PSF data
% 
% % Save the data
% ieSaveSIDataFile(psf,wave,umPerSample,'SI-pillBox');
