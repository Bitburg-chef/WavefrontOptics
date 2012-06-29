% v_wvfAstigmatism
%
% Compute the PSF for various astigmatism and blur levels.  This
% illustrates the effect of Zernicke coefficients 4 and 5
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

%% Calculate point spread for wavelength defocus

%% Set up default parameters structure with diffraction limited default
wvfP = wvfCreate;
wvfParams = wvfComputePSF(wvfP);
z = wvfGet(wvfParams,'zcoeffs');
z4 = -0.5:0.5:0.5; z5 = -0.5:0.5:0.5;
[Z4,Z5] = meshgrid(z4,z5);
Zvals = [Z4(:), Z5(:)];

%% Make the plot
vcNewGraphWin;
for ii=1:size(Zvals,1)
    subplot(3,3,ii)
    z(4) = Zvals(ii,1); z(5) = Zvals(ii,2);
    wvfParams = wvfSet(wvfParams,'zcoeffs',z);
    wvfParams = wvfComputePSF(wvfParams);
    wvfPlot(wvfParams,'2d psf space','um',1,maxUM);
    title(sprintf('W4 = %.1f W5 == %.1f\n',z(4),z(5)));
end





