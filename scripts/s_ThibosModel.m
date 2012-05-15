%% Explore the Thibos data and model for human point spread functions
%
% The Zernicke coefficients for the human eye at different pupil sizes have
% been estimated and the data distributed by Thibos and the gang in
% Indiana.  Brainard snagged them along with the Virtual Eyes demo stuff.
% We re-wrote and integrated with the toolbox here.
%
% This routine loads the data and plots them using our methods.
%
% See also:  VirtualEyesDemo script therein.
%
% Copyright Wavefront Toolbox Team, 2012

%% Initialize
s_initISET

maxUM = 40;

%%  Load the statistical wavefront properties for a 6.0mm pupil diameter

% The Zernike coefficients describing the wavefront aberrations are each
% distributed as a Gaussian. There is some covariance between these
% coefficients.  The covariance is summarized in the variable S.  The mean
% values across a large sample of eyes measured by Thibos and gang are in
% the variable sample_mean.

pupilMM = 4.5;
[sample_mean S] = vwfLoadHuman(pupilMM);


%% Plot the means and covariance (not)

% Put some utilities in vcNewGraphWin for different screen formats.  Use
% this information, say for a 3,1 on the left side, or a 1,3 towards the
% top.

ss = get(0,'ScreenSize');
set(vcNewGraphWin,'position',[ 76   152   369   757]);

subplot(3,1,1)
plot(sample_mean,'--o'); grid on
xlabel('Zernike polynomial number')
ylabel('Coefficient value')

subplot(3,1,2)
imagesc(S);
axis image, title('Covariance')
colormap(hot); colorbar

N = 10;

% Calculate sample eyes using the multivariate normal distribution
% Each column of Zcoeffs is an example person
Zcoeffs = ieMvnrnd(sample_mean,S,N)';  % each row of R is a vector of Zernike coeffs

% Show 'em
subplot(3,1,3)
plot(Zcoeffs); grid on
xlabel('Zernike polynomial number')
ylabel('Coefficient value')

%% Examine a single PSF for a single example

% Allocate space and fill in the lower order coefficients
z = zeros(65,1);
whichSample = 3;
z(1:13) = Zcoeffs(1:13,whichSample);

% Create the typical subject
thisGuy = wvfCreate;                                 % Initialize
thisGuy = wvfSet(thisGuy,'zcoeffs',z);               % Zernike
thisGuy = wvfSet(thisGuy,'measuredpupil',pupilMM);         % Data
thisGuy = wvfSet(thisGuy,'calculatedpupil',pupilMM); % What we calculate
thisGuy = wvfSet(thisGuy,'nominalfocuswl',550);
thisGuy = wvfSet(thisGuy,'wavelength',[500 50 3]);  % SToWls format
thisGuy = wvfComputePSF(thisGuy);

%% Plot the PSFs for the distinct wavelengths
wave = wvfGet(thisGuy,'wave');
nWave= wvfGet(thisGuy,'nwave');
% for ii=1:nWave
%     vcNewGraphWin; 
%     mesh(thisGuy.psf{ii}); title(sprintf('%d nm',wave(ii)));
% end
vcNewGraphWin([],'tall');
for ii=1:nWave
    subplot(nWave,1,ii)
    wvfPlot(thisGuy,'image psf space','um',ii,maxUM);
    title(sprintf('%d nm',wave(ii)));
end

% close all

%% Not currently used, but could be interesting

% vcNewGraphWin;
% wvfPlot(thisGuy,'2d pupil phase space','mm',pupilfuncrangeMM);


%% Calculate the PSFs from the coeffcients

% Allocate space
z = zeros(65,1);

f = vcNewGraphWin([],'tall');
% For some examples
pupilfuncrangeMM = 6;
whichSubjects = 1:3:N;
nSubjects = length(whichSubjects);
for ii = 1:nSubjects
    
    % whichSubject = 1;
    z(1:13) = Zcoeffs(1:13,whichSubjects(ii));
    
    % Place them
    thisGuy = wvfCreate;                                 % Initialize
    thisGuy = wvfSet(thisGuy,'zcoeffs',z);  % Zernike
    thisGuy = wvfSet(thisGuy,'measuredpupil',pupilMM);  % Data
    thisGuy = wvfSet(thisGuy,'calculatedpupil',pupilMM);% What we calculate
    thisGuy = wvfComputePSF(thisGuy);

    subplot(nSubjects,1,ii)
    wvfPlot(thisGuy,'image psf space','um',1,maxUM);
    title(sprintf('Subject %d\n',ii))
end


