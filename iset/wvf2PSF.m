function [psf, wave, umPerSample] = wvf2PSF(wvfP)
% Convert a wvf structure to an ISET shift-invariant PSF
%
%    [psf, wave, umPerSample] = wvf2PSF(wvfP)
%
% For each wavelength in wvf, compute the PSF with proper units and
% place it in an ISET shift-invariant PSF format that can be used for human
% optics simulation.
%
% Example:
%    pupilMM = 6; wls = [400 50 7]; wvfP = wvfLoadHuman(pupilMM,wave);
%    [psf, wave, umPerSample] = wvf2PSF(wvfP)
%    fName = sprintf('psfSI-%s',wvfGet(wvfP,'name'));
%    ieSaveSIDataFile(psf,wave,umPerSample,fName);

% Copyright Imageval 2012

%% Parameters
if ieNotDefined('wvfP'), error('wvf parameters required.'); end
wave = wvfGet(wvfP,'wave');
nWave = wvfGet(wvfP,'nwave');

% PSF - number of pixels and spacing in microns between samples
nPix = 256;                 % Number of pixels
umPerSample = [0.25,0.25];  % In microns

wvfP = wvfComputePSF(wvfP);
samp = wvfGet(wvfP,'samples space','um');
samp = samp(:);

iSamp = (1:nPix)*umPerSample(1);
iSamp = iSamp - mean(iSamp);
iSamp = iSamp(:);

psf = zeros(nPix,nPix,nWave);

% I am worried that we don't have the sampling right for the wave
% PSF data
% vcNewGraphWin([],'upperleft'); unit = 'um';  pRange = 40;
for ii=1:nWave,
    thisPSF = wvfGet(wvfP,'psf',ii);
    psf(:,:,ii) = interp2(samp,samp',thisPSF,iSamp,iSamp');
    % wvfPlot(wvfP,'image psf space',unit,ii,pRange)
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
