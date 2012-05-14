function siPSF = wvf2PSF(wvfP)
% Convert a wvf structure to an ISET shift-invariant PSF
%
%    siPSF = wvf2PSF(wvfP)
%
% For each wavelength in wvf, compute the PSF with proper units and
% place it in an ISET shift-invariant PSF format that can be used for human
% optics simulation.
%
% Copyright Imageval 2012

if ieNotDefined('wvfP'), error('wvf parameters required.'); end


psf = wvfComputePSF(wvfP);

pupilMM = 4.5;
load('IASstats45','sample_mean');
zCoeffs = zeros(65,1);
zCoeffs(1:length(sample_mean)) = sample_mean;

wvfP = wvfCreate;
wvfP = wvfSet(wvfP,'zcoeffs',zCoeffs);               % Zernike
wvfP = wvfSet(wvfP,'measuredpupil',6);         % Data
wvfP = wvfSet(wvfP,'calculatedpupil',pupilMM); % What we calculate
wvfP = wvfSet(wvfP,'nominalfocuswl',550);
wvfP = wvfSet(wvfP,'wavelength',[400 50 7]);  % SToWls format
wvfP = wvfComputePSF(wvfP);

end


%% From ISET script on how to create an SI data file

% Now, write out a file containing the relevant point spread function
% data, along with related variables.
umPerSample = [0.25,0.25];                % Sample spacing

% Point spread is a little square in the middle of the image
h = zeros(128,128); h(48:79,48:79) = 1; h = h/sum(h(:));
for ii=1:length(wave), psf(:,:,ii) = h; end     % PSF data

% Save the data
ieSaveSIDataFile(psf,wave,umPerSample,'SI-pillBox');
