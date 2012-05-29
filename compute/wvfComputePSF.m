function wvfP = wvfComputePSF(wvfP)
% Compute the monochromatic psf over the calculated pupil
%
%    wvfP = wvfComputePSF(wvfP)
%
% wvfP: A structure of wavefront parameters. See wvfCreate for the default
%       fields
%
% The point spread function is computed for each of the wavelengths listed
% in the input wvfP structure. The PSF computation is based on 10 orders of
% Zernike coefficients specified to the OSA standard.
%
% The computation includes the Stiles-Crawford Effect (SCE) if specified by
% the wvfP structure.
%
% The psfs at each wavelength are scaled so that each sums to unity.  If
% SCE is specified, you can multiply by sceFrac at that wavelength to
% incorporate the fraction of light absorbed.
%
% Required input fields for wvfP struct -
%  (see comment in wvfComputePupilFunction for more details)
%
%   zcoeffs -           Zernike coefficients.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%   caclpupilsize -     Size over which returned pupil function is calculated, in MM.
%   wls -               Column vector of wavelengths over which to compute, in NANOMETERS.
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   fieldSampleSizeMMperPixel - Size in mm of each pixel of the pupil
%                       function.
%   sizeOfFieldMM -     Size of square image over which the pupile function is computed in MM.
%
% Optional input fields for wvfP struct
%   sceParams -         Parameter structure for Stiles-Crawford correction.  If missing or set to empty,
%                       no correction and is set to empty on return.
%
% Output fields set in wvfP struct
%   psf -               Calcuated polychromatic psf. Third dimension of returned matrix indexes wavelength.
%   pupilfunc -         Calculated pupil function.  Third dimension of returned matrix indexes wavelength
%   arcminperpix -      Arc minutes per pixel for returned psfs.
%   strehl -            Strehl ratio of psf at each wavelength.  If SCE correction is specified, the returned
%                       strehl ratio is to the diffraction limited psf with the same SCE assumed.
%   sceFrac -           Fraction of light actually absorbed when SCE is taken into account, at each wavelength.
%   areapix -           Number of pixels within the computed pupil aperture at each wavelength
%   areapixapod -       Number of pixels within the computed pupil aperture at each wavelength,
%                       multiplied by the Stiles-Crawford aopdization.
%   defocusMicrons -    Defocus added in to zcoeffs(4) at each wavelength, in microns.
%
% Note: These calculations only account for longitudinal chromatic
% aberration (LCA); they do NOT incoporate transverse chromatic aberration
% (TCA) or off-axis at all.
%
% NOTES
% 1. When zcoeffs is a vector of zeros, the output monochromatic and
% polychromatic PSFs are diffraction-limited. and longitudinal chromatic
% aberration (and the SCE if specified).
%
% 2. There are three ways that defocus can be specified to this routine.
%   A. First, the value of the 4th zernike coefficient is included.
%   B. Second, value is computed based on the differencebetween nominalFocusWl
%   and the wavelength for which the PSF is being computed.
%   C. Third, the value of defocusDiopters is added in directly.
%
%  These are redundant, it is often most convenient to think in terms of
%  one of the three ways.
%  KP 3/12/12: Only A and C are redundant? B is the result of longitudinal
%  chromatic aberration (LCA) for using non-nominal wavelength to measure
%  PSF. A and C are defocus due to monochromatic aberrations.
%
% See also: wvfComputeConePSF, wvfComputePupilFunction, sceGetParamsParams,
% wvfGetDefocusFromWavelengthDifference
%
% Based on code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.
% 9/5/11  dhb      Rename.  Rewrite for wvfP i/o.
%
% Copyright Wavefront Toolbox Team 2012

%% Check optional fields
if (~isfield(wvfP,'sceParams') || isempty(wvfP.sceParams))
    wvfP.sceParams = [];
end


%% Make sure pupil function is computed
wvfP = wvfComputePupilFunction(wvfP);

wave = wvfGet(wvfP,'wave');
nWave = wvfGet(wvfP,'nwave');

psf = cell(nWave,1);
areapix = zeros(nWave,1);
areapixapod = zeros(nWave,1);
strehl  = zeros(nWave,1);
sceFrac = zeros(nWave,1);

% We will set tmpWvfParams to a single wavelength and loop.
w = waitbar(0,'Pupil functions');
for wl = 1:nWave
    
    waitbar(wl/nWave,w,sprintf('Psf %d',wave(wl)));    
    
    % Convert to the pupil function to the PSF  Just this simple.
    % Scale so that psf sums to unity.
    pupilfunc{wl} = wvfGet(wvfP,'pupil function',wl);
    amp = fft2(pupilfunc{wl});
    inten = (amp .* conj(amp));   %intensity
    psf{wl} = real(fftshift(inten));
    psf{wl} = psf{wl}/sum(sum(psf{wl}));
    
    % Get strehl ratio, based on Heidi's code
    % strehl(wl) = max(max(psf{wl}))./(areapixapod(wl)^2);
    
    % The variable sceFrac tells you how much light
    % is lost if you correct for the SCE.
    % sceFrac(wl) = (areapixapod(wl)/areapix(wl));
end
close(w);  %Waitbar

%% Set output fields - these are vectors or 3D matrices
wvfP = wvfSet(wvfP,'psf',psf);
% wvfP = wvfSet(wvfP,'strehl',strehl);
% wvfP = wvfSet(wvfP,'areapix',areapix);
% wvfP = wvfSet(wvfP,'areapixapod',areapixapod);

% Shouldn't we be setting areapixapod and areapix, too?

% Derived, not needed? wvfP.arcminperpix = arcminperpix;
% wvfP.strehl = strehl;
% wvfP.sceFrac = sceFrac;
% Don't understand yet. (BW)
% wvfP.areapix = areapix;
% wvfP.areapixapod = areapixapod;

end


