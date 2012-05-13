%% Explore the Thibos data and model for human point spread functions
%
% ieMvnrnd call
% Thibos data correlation matrix.
%
% See the VirtualEyesDemo script therein.
%

%%
s_initISET

%%  Load the statistical wavefront properties for a 6.0mm pupil diameter

% The Zernike coefficients describing the wavefront aberrations are each
% distributed as a Gaussian. There is some covariance between these
% coefficients.  The covariance is summarized in the variable S.  The mean
% values across a large sample of eyes measured by Thibos and gang are in
% the variable sample_mean.

% These are the summary statistics for a 6.0mm pupil
load('IASstats60','S','sample_mean');

% Plot the means and covariance (not)
vcNewGraphWin;
subplot(1,2,1)
plot(sample_mean,'--o'); grid on
xlabel('Zernike polynomial number')
ylabel('Coefficient value')

subplot(1,2,2)
imagesc(S);
axis image, title('Covariance')
colormap(hot); colorbar

%% Calculate sample eyes using the multivariate normal distribution

% Each column of Zcoeffs is an example person
N = 10;
Zcoeffs = ieMvnrnd(sample_mean,S,N)';  % each row of R is a vector of Zernike coeffs

% Show 'em
vcNewGraphWin;
plot(Zcoeffs); grid on
xlabel('Zernike polynomial number')
ylabel('Coefficient value')

%% Examine a single PSF for a single example

% Allocate space and fill in the lower order coefficients
z = zeros(65,1);
z(1:13) = Zcoeffs(1:13,1);

% Create the typical subject
thisGuy = wvfCreate;                                 % Initialize
thisGuy = wvfSet(thisGuy,'zcoeffs',z);  % Zernike
thisGuy = wvfSet(thisGuy,'measuredpupil',6);  % Data
thisGuy = wvfSet(thisGuy,'calculatedpupil',6);% What we calculate
thisGuy = wvfSet(thisGuy,'nominalfocuswl',550);

% Not sure what this is used for
% pupilfuncrangeMM = 6;
psffuncmaxMin = 1/3;

% Compute for the current wavelenth
thisGuy550 = wvfSet(thisGuy,'wavelength',550);
thisGuy550 = wvfComputePSF(thisGuy550);

vcNewGraphWin;
wvfPlot(thisGuy550,'2d psf angle normalized','deg',psffuncmaxMin);

%% Now for 450
thisGuy450 = wvfSet(thisGuy,'wavelength',450);
thisGuy450 = wvfComputePSF(thisGuy450);
vcNewGraphWin;
wvfPlot(thisGuy450,'2d psf angle normalized','deg',psffuncmaxMin);

%% And now for 650
thisGuy650 = wvfSet(thisGuy,'wavelength',650);
thisGuy650 = wvfComputePSF(thisGuy650);
vcNewGraphWin;
wvfPlot(thisGuy650,'2d psf angle normalized','deg',psffuncmaxMin);

%% Not currently used, but could be interesting

% vcNewGraphWin;
% wvfPlot(thisGuy,'2d pupil phase space','mm',pupilfuncrangeMM);


%% Calculate the PSFs from the coeffcients

% Allocate space
z = zeros(65,1);

% For some examples
for whichSubject = 1:2:N
    % whichSubject = 1;
    z(1:13) = Zcoeffs(1:13,whichSubject);
    
    % Place them
    thisGuy = wvfCreate;                                 % Initialize
    thisGuy = wvfSet(thisGuy,'zcoeffs',z);  % Zernike
    thisGuy = wvfSet(thisGuy,'measuredpupil',6);  % Data
    thisGuy = wvfSet(thisGuy,'calculatedpupil',6);% What we calculate
    
    thisGuy = wvfComputePSF(thisGuy);
    pupilfuncrangeMM = 6;
    psffuncmaxMin = 1/3;
    
    % vcNewGraphWin;
    % wvfPlot(thisGuy,'2d pupil phase space','mm',pupilfuncrangeMM);
    
    vcNewGraphWin;
    wvfPlot(thisGuy,'2d psf angle normalized','deg',psffuncmaxMin);
    title(sprintf('Subject %d\n',whichSubject))
end


