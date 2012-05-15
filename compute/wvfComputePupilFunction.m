function [wvfP, phase, A] = wvfComputePupilFunction(wvfP, varargin)
% Compute the monochromatic pupil fuction
%
%    [wvfP,phase,amplitude] = wvfComputePupilFunction(wvfP, [idxWave=1])
%
% The pupil function is a complex number that represents the amplitude and
% phase of the wavefront across the pupil.  The returned pupil function at
% a specific wavelength (in microns) is
%
%    pupilF = A exp(-1i 2 pi (phase/wavelength));
%
% The amplitude is calculated entirely based on the assumed properties of
% the Stiles-Crawford effect.
%
% The pupil function is related to the PSF by the Fourier transform. See J.
% Goodman, Intro to Fourier Optics, 3rd ed, p. 131. (MDL)
%
% These functions are calculated for 10 orders of Zernike coeffcients specified to
% the OSA standard. Includes SCE (Stiles-Crawford Effect) if specified.
% The SCE is modeled as an apodization filter (a spatially-varying amplitude
% attenuation) to the pupil function. In this case, it is a decaying exponential.
% (MDL)
%
% Required input fields for wvfP struct
%   zcoeffs -           Zernike coefficients. Expects 65 coefficients numbered with the osa j index.
%                       These are up to the 10th order.  Coefficients zcoeffs(1) and zcoeffs(2) are tip and
%                       tilt, and are not entered into the calculations (i.e. treated as zero ). zcoeffs(3)
%                       through zcoeffs(5) are astigmatism and defocus.  You can pass fewer than the full 65
%                       coefficients, in which case the trailing coefficients are set to zero.
%
%   measpupilMM -       Size of pupil characterized by the coefficients, in
%     MM. This is how large the physical pupil is and determines the
%     normalized scaling of rho and the Zernike polynomals.
%
%   caclpupilMM -       Size over which returned pupil function is
%     calculated, in MM. Must be smaller than measpupilMM. The pupil function
%     is set to zero outside this radius.
%
%   wls -               Wavelength to compute for, in NM.  Can only pass one wavelenth, despite plural in the name.
%                       This is because wvfComputePupilFunction(tmpwvfParams) is passed a temporary wvf by wvfComputePSF
%                       which only has 1 wavelength. wvfComputePSF handles the loop through the vector of wls
%                       that the original wvfP contains.
%
%   fieldSampleSizeMMperPixel - Size in mm of each pixel of the pupil
%                       function.
%
%   sizeOfFieldMM -     Size of square image over which the pupil function is computed in MM.
%                       Setting this larger than the calculated pupil size prevents undersampling
%                       of the PSF that will ultimately be computed from the pupil function.
%                       MDL: I guess this is here so that there is some
%                       zero-padding of the calculated pupil function. This
%                       allows the calculated PSF (Fouier transform of the
%                       pupil) to be sampled at a higher spatial
%                       resolution.
%
% Optional input fields for wvfP struct
%   sceParams -         Parameter structure for Stiles-Crawford correction.  If missing or set to empty,
%                       no correction and is set to empty on return.  See sceGetParamsParams.
%
% Output fields set in wvfP struct
%   pupilfunc -     Calcuated pupil function
%
% PROGRAMMING NOTE:  The notion of pixel isn't so good.  We need to replace
% it with a measure that has a clear physical description throughout.  If
% it is always sample, then the sample must have a spatial size in um or
% mm. The good news is I think this is the last item that doesn't have an
% easily identified physical unit.
%
% 3/16/12 MDL has changed the sizeOfFieldPixels to
% fieldSampleSizeMMPerPixel as the complement to sizeOfFieldMM. This should
% reinforce the notion of a physical unit behind the code, where the size
% of the field can be specified in MM, as well as a sampling rate of
% MM/pixel. The code then computes the pixel value from that.
%
% All aberrations other than defocus (including astigmatism) are assumed to
% be constant with wavelength, as variation with wavelength in other
% aberrations is known to be small.
%
% Transverse chromatic aberration (TCA), which is a wavelength dependent tip
% or tilt, has also not been included.
%
% Dividing the psf computed from the returned pupil function by areapix
% (or areapixapod) squared effects a normalization so that the peak is
% the strehl ratio.
%
% See also: wvfComputePSF, sceGetParamsParams.
%
% Code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine.
% Reformat comments. 
% 9/5/11  dhb      Rewrite for wvfP struct i/o.  Rename.
%
% (c) Wavefront Toolbox Team 2011, 2012

%% Parameter checking
if ieNotDefined('wvfP'), error('wvfP required'); end

% Check pupil size issue
calcPupilSizeMM = wvfGet(wvfP,'calculated pupil','mm');
measPupilSizeMM = wvfGet(wvfP,'measured pupil','mm');
if (calcPupilSizeMM > measPupilSizeMM)
    error('Calculation pupil (%.2f mm) must not exceed measurement pupil (%.2f mm).', ...
        calcPupilSizeMM, measPupilSizeMM);
end

% Handle case where not all 65 coefficients are passed
c = zeros(65,1);
c(1:length(wvfGet(wvfP,'zcoeffs'))) = wvfGet(wvfP,'zcoeffs');

%% Calculate pupil function for each wavelength

% Convert wavelengths in nanometers to wavelengths in microns
% wlInUM = wvfP.wls/1000;
% Sanity check
waveUM = wvfGet(wvfP,'wave','um');
waveNM = wvfGet(wvfP,'wave','nm');
nWave = wvfGet(wvfP,'n wave');

for ii=1:nWave
    thisWave = waveNM(ii);
    
    % Set SCE correction params, if desired
    xo  = wvfGet(wvfP,'scex0');
    yo  = wvfGet(wvfP,'scey0');
    rho = wvfGet(wvfP,'sce rho');
    
    % set up pupil coordinates to compute A and phase
    % Create arrays that represent the length coordinates of the pupil
    nPixels = wvfGet(wvfP,'npixels');
    fieldSizeMM = wvfGet(wvfP,'fieldsize','mm',ii);
    pupilPos = (0:(nPixels-1))*(fieldSizeMM/nPixels)-fieldSizeMM/2;
    [xpos ypos] = meshgrid(pupilPos);
    
    % Set up the amplitude of the pupil function.
    % This appears to depend entirely on the SCE correction
    if all(rho) == 0, A=ones(nPixels);
    else
        % Get the wavelength-specific value of rho for the Stiles-Crawford
        % effect.
        rho      = wvfGet(wvfP,'sce rho',thisWave);
        
        % For the x,y positions within the pupil, the value of rho is used to
        % set the amplitude.  I guess this is where the SCE stuff matters.  We
        % should have a way to expose this for teaching and in the code.
        
        % 3/9/2012, MDL: Removed nested for loop for calculating the
        % SCE. Note previous code had x as rows of matrix, y as columns of
        % matrix. This has been corrected so that x is columns, y is rows.
        A=10.^(-rho*((xpos-xo).^2+(ypos-yo).^2));
    end
    
    % 3/9/2012, MDL: Removed nested for loop for calculating the
    % phase. Note previous code had x as rows of matrix, y as columns of
    % matrix. This has been corrected so that x is columns, y is rows.
    % Also renamed "angle" to "theta" since angle is a built-in MATLAB function.
    % Also redefined angle in terms of natural polar coordinates instead of
    % prevous definition:
    %        if (xpos==0 && ypos>0),     angle = 3.1416/2;
    %        elseif(xpos==0 && ypos<0),  angle = -3.1416/2;
    %        elseif(xpos==0 && ypos==0), angle = 0;
    %        elseif(xpos>0),             angle = atan(ypos/xpos);
    %        else                        angle= 3.1416 + atan(ypos/xpos);
    %        end
    norm_radius = (sqrt(xpos.^2+ypos.^2))/(measPupilSizeMM/2);
    theta = atan2(ypos,xpos);
    phase = 0 + ...
        c(5) .* sqrt(6).*norm_radius.^2 .* cos(2 .* theta) + ...
        c(3) .* sqrt(6).*norm_radius.^2 .* sin(2 .* theta) + ...
        c(4) .* sqrt(3).*(2 .* norm_radius.^2 - 1) + ...
        c(9) .*sqrt(8).* norm_radius.^3 .* cos(3 .* theta) + ...
        c(6) .*sqrt(8).* norm_radius.^3 .* sin(3 .* theta) + ...
        c(8) .*sqrt(8).* (3 .* norm_radius.^3 - 2 .* norm_radius) .* cos(theta) + ...
        c(7) .*sqrt(8).* (3 .* norm_radius.^3 - 2 .* norm_radius) .* sin(theta) + ...
        c(14) .* sqrt(10).*norm_radius.^4 .* cos(4 .* theta) + ...
        c(10) .* sqrt(10).*norm_radius.^4 .* sin(4 .* theta) + ...
        c(13) .* sqrt(10).*(4 .* norm_radius.^4 - 3 .* norm_radius.^2) .* cos(2 .* theta) + ...
        c(11) .* sqrt(10).*(4 .* norm_radius.^4 - 3 .* norm_radius.^2) .* sin(2 .* theta) + ...
        c(12) .* sqrt(5).*(6 .* norm_radius.^4 - 6 .* norm_radius.^2 + 1)+...
        c(20) .* 2.*sqrt(3).*norm_radius.^5 .* cos(5 .* theta) + ...
        c(15) .*2.*sqrt(3).* norm_radius.^5 .* sin(5 .* theta) + ...
        c(19) .* 2.*sqrt(3).*(5 .* norm_radius.^5 - 4 .* norm_radius.^3) .* cos(3 .* theta) + ...
        c(16) .*2.*sqrt(3).* (5 .* norm_radius.^5 - 4 .* norm_radius.^3) .* sin(3 .* theta) + ...
        c(18) .*2.*sqrt(3).* (10 .* norm_radius.^5 - 12 .* norm_radius.^3 + 3 .* norm_radius) .* cos(theta) + ...
        c(17) .*2.*sqrt(3).* (10 .* norm_radius.^5 - 12 .* norm_radius.^3 + 3 .* norm_radius) .* sin(theta) + ...
        c(27) .*sqrt(14).* norm_radius.^6 .* cos(6 .* theta) + ...
        c(21) .*sqrt(14).*norm_radius.^6 .* sin(6 .* theta) + ...
        c(26) .*sqrt(14).*(6 .* norm_radius.^6 - 5 .* norm_radius.^4) .* cos(4 .* theta) + ...
        c(22) .*sqrt(14).*(6 .* norm_radius.^6 - 5 .* norm_radius.^4) .* sin(4 .* theta) + ...
        c(25) .*sqrt(14).* (15 .* norm_radius.^6 - 20 .* norm_radius.^4 + 6 .* norm_radius.^2) .* cos(2 .* theta) + ...
        c(23) .*sqrt(14).*(15 .* norm_radius.^6 - 20 .* norm_radius.^4 + 6 .* norm_radius.^2) .* sin(2 .* theta) + ...
        c(24) .*sqrt(7).* (20 .* norm_radius.^6 - 30 .* norm_radius.^4 + 12 .* norm_radius.^2 - 1)+...
        c(35) .*4.* norm_radius.^7 .* cos(7 .* theta) + ...
        c(28) .*4.* norm_radius.^7 .* sin(7 .* theta) + ...
        c(34) .*4.* (7 .* norm_radius.^7 - 6 .* norm_radius.^5) .* cos(5 .* theta) + ...
        c(29) .*4.* (7 .* norm_radius.^7 - 6 .* norm_radius.^5) .* sin(5 .* theta) + ...
        c(33) .*4.* (21 .* norm_radius.^7 - 30 .* norm_radius.^5 + 10 .* norm_radius.^3) .* cos(3 .* theta) + ...
        c(30) .*4.* (21 .* norm_radius.^7 - 30 .* norm_radius.^5 + 10 .* norm_radius.^3) .* sin(3 .* theta) + ...
        c(32) .*4.* (35 .* norm_radius.^7 - 60 .* norm_radius.^5 + 30 .* norm_radius.^3 - 4 .* norm_radius) .* cos(theta) + ...
        c(31) .*4.* (35 .* norm_radius.^7 - 60 .* norm_radius.^5 + 30 .* norm_radius.^3 - 4 .* norm_radius) .* sin(theta) +...
        c(44) .*sqrt(18).* norm_radius.^8 .* cos(8 .* theta) + ...
        c(36) .*sqrt(18).* norm_radius.^8 .* sin(8 .* theta) + ...
        c(43) .*sqrt(18).* (8 .* norm_radius.^8 - 7 .* norm_radius.^6) .* cos(6 .* theta) + ...
        c(37) .*sqrt(18).* (8 .* norm_radius.^8 - 7 .* norm_radius.^6) .* sin(6 .* theta) + ...
        c(42) .*sqrt(18).* (28 .* norm_radius.^8 - 42 .* norm_radius.^6 + 15 .* norm_radius.^4) .* cos(4 .* theta) + ...
        c(38) .*sqrt(18).* (28 .* norm_radius.^8 - 42 .* norm_radius.^6 + 15 .* norm_radius.^4) .* sin(4 .* theta) + ...
        c(41) .*sqrt(18).* (56 .* norm_radius.^8 - 105 .* norm_radius.^6 + 60 .* norm_radius.^4 - 10 .* norm_radius.^2) .* cos(2 .* theta) + ...
        c(39) .*sqrt(18).* (56 .* norm_radius.^8 - 105 .* norm_radius.^6 + 60 .* norm_radius.^4 - 10 .* norm_radius.^2) .* sin(2 .* theta) + ...
        c(40) .*3.* (70 .* norm_radius.^8 - 140 .* norm_radius.^6 + 90 .* norm_radius.^4 - 20 .* norm_radius.^2 + 1) + ...
        c(54) .*sqrt(20).* norm_radius.^9 .* cos(9 .* theta) + ...
        c(45) .*sqrt(20).* norm_radius.^9 .* sin(9 .* theta) + ...
        c(53) .*sqrt(20).* (9 .* norm_radius.^9 - 8 .* norm_radius.^7) .* cos(7 .* theta) + ...
        c(46) .*sqrt(20).* (9 .* norm_radius.^9 - 8 .* norm_radius.^7) .* sin(7 .* theta) + ...
        c(52) .*sqrt(20).* (36 .* norm_radius.^9 - 56 .* norm_radius.^7 + 21 .* norm_radius.^5) .* cos(5 .* theta) + ...
        c(47) .*sqrt(20).* (36 .* norm_radius.^9 - 56 .* norm_radius.^7 + 21 .* norm_radius.^5) .* sin(5 .* theta) + ...
        c(51) .*sqrt(20).* (84 .* norm_radius.^9 - 168 .* norm_radius.^7 + 105 .* norm_radius.^5 - 20 .* norm_radius.^3) .* cos(3 .* theta) + ...
        c(48) .*sqrt(20).* (84 .* norm_radius.^9 - 168 .* norm_radius.^7 + 105 .* norm_radius.^5 - 20 .* norm_radius.^3) .* sin(3 .* theta) + ...
        c(50) .*sqrt(20).* (126 .* norm_radius.^9 - 280 .* norm_radius.^7 + 210 .* norm_radius.^5 - 60 .* norm_radius.^3 + 5 .* norm_radius) .* cos(theta) + ...
        c(49) .*sqrt(20).* (126 .* norm_radius.^9 - 280 .* norm_radius.^7 + 210 .* norm_radius.^5 - 60 .* norm_radius.^3 + 5 .* norm_radius) .* sin(theta) + ...
        c(65) .*sqrt(22).* norm_radius.^10 .* cos(10 .* theta) + ...
        c(55) .*sqrt(22).* norm_radius.^10 .* sin(10 .* theta) + ...
        c(64) .*sqrt(22).* (10 .* norm_radius.^10 - 9 .* norm_radius.^8) .* cos(8 .* theta) + ...
        c(56) .*sqrt(22).* (10 .* norm_radius.^10 - 9 .* norm_radius.^8) .* sin(8 .* theta) + ...
        c(63) .*sqrt(22).* (45 .* norm_radius.^10 - 72 .* norm_radius.^8 + 28 .* norm_radius.^6) .* cos(6 .* theta) + ...
        c(57) .*sqrt(22).* (45 .* norm_radius.^10 - 72 .* norm_radius.^8 + 28 .* norm_radius.^6) .* sin(6 .* theta) + ...
        c(62) .*sqrt(22).* (120 .* norm_radius.^10 - 252 .* norm_radius.^8 + 168 .* norm_radius.^6 - 35 .* norm_radius.^4) .* cos(4 .* theta) + ...
        c(58) .*sqrt(22).* (120 .* norm_radius.^10 - 252 .* norm_radius.^8 + 168 .* norm_radius.^6 - 35 .* norm_radius.^4) .* sin(4 .* theta) + ...
        c(61) .*sqrt(22).* (210 .* norm_radius.^10 - 504 .* norm_radius.^8 + 420 .* norm_radius.^6 - 140 .* norm_radius.^4 + 15 .* norm_radius.^2) .* cos(2 .* theta) + ...
        c(59) .*sqrt(22).* (210 .* norm_radius.^10 - 504 .* norm_radius.^8 + 420 .* norm_radius.^6 - 140 .* norm_radius.^4 + 15 .* norm_radius.^2) .* sin(2 .* theta) + ...
        c(60) .*sqrt(11).* (252 .* norm_radius.^10 - 630 .* norm_radius.^8 + 560 .* norm_radius.^6 - 210 .* norm_radius.^4 + 30 .* norm_radius.^2 - 1);
    
    % Here is the pupil function
    pupilfunc = A.*exp(-1i * 2 * pi * phase/waveUM(ii));
    
    % Set values outside the radius to to 0
    pupilfunc(norm_radius > calcPupilSizeMM/measPupilSizeMM)=0;
    
    % Attach the function the the proper wavelength slot
    wvfP = wvfSet(wvfP,'pupilfunc',pupilfunc,ii);
    
end

end

