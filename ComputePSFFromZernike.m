function [psf,arcminperpix,strehl,sceFrac,areapix,areapixapod] = ComputePSFFromZernike(zcoeffs,measpupilMM,calcpupilMM,wls,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,sceParams)
% [psf,arcminperpix,strehl,sceFrac.areapix,areapixapod] = ComputePSFFromZernike(zcoeffs,measpupilMM,calcpupilMM,wls,nominalFocusWl,defocusDiopters,sizeOfFieldPixels,sizeOfFieldMM,[sceParams])
%
% Computes the monochromatic psf over the calculated pupil size for 10 orders of Zernike
% coeffcients specified to the OSA standard. Includes SCE (Stiles-Crawford Effect) if specified.
%
% The psfs at each wavelength are scaled so that each sums to unity.  If SCE is specified,
% you can multiply by sceFrac at that wavelength to get the fraction of light absorbed.
%
% Inputs - see comment in ComputePupilFunctionFromZernike for more details.
%   zcoeffs -           Zernike coefficients.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%   caclpupilsize -     Size over which returned pupil function is calculated, in MM.
%   wls -               Column vector of wavelengths over which to compute, in NANOMETERS.
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   sizeOfFieldPixels - Linear size of square image over which the pupil function is computed.
%   sizeOfFieldMM -     Size of square image over which the pupile function is computed in MM.
%   sceParams -         Parameter structure for Stiles-Crawford correction. Default, no correction.
%
% Outputs
%   psf -               Calcuated polychromatic psf. Third dimension of returned matrix indexes wavelength.                      
%   arcminperpix -      Arc minutes per pixel for returned psfs.
%   strehl -            Strehl ratio of psf at each wavelength.
%   sceFrac -           Fraction of light actually absorbed when SCE is taken into account, at each wavelength.
%   areapix -           Number of pixels within the computed pupil aperture at each wavelength
%   areapixapod -       Number of pixels within the computed pupil aperture at each wavelength,
%                       multiplied by the Stiles-Crawford aopdization.
%
% Note: If SCE is included the returned polychromatic psfs will generally integrate to a number
% smaller than 1, this reflects light assumed not to have been caught and guided by the receptors.
%
% Note: These calculations only account for lognitudinal chromatic aberration (LCA), and do not
% incoporate transverse chromatic aberration (TCA).
%
% Note: If this function is called with a zcoeffs vector that is composed of only zeros 
% the output monochromatic and polychromatic PSFs will be
% limited only by diffraction and longitudinal chromatic aberration (and the SCE if specified). 
%
% See also: ComputePupilFunctionFromZernike, GetStilesCrawfordParams.
%
% Code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.

%% Check optional args
if (nargin < 9 || isempty(sceParams))
    sceParams = [];
end

%% Handle defocus relative to reference wavelength.
[microns,diopters] = GetDefocusFromWavelengthDifference(wls,nominalFocusWl,defocusDiopters,measpupilMM);

% Store the original value of defocus (zcoeffs(4)) since it will get changed
% later when changing wavelength
doriginal = zcoeffs(4);

%% Get the pupil function, correcting defocus at each wavelength
%
% The value of arcminperpix is set depending on a wavelength chosen
% to set the scale.  Here is what Heidi says about this process,
% and about the renormalization of scale for each call to
% ComputePupilFunctionFromZernike:
%
%   The setting of the scale and then the change in the
%   loop with wavelength is to make sure that the scale on the psfs stays
%   the same as the wavelength is changed - since the scale is related to
%   lamda/(linear dimension of matrix) - this is just related to the way
%   the dimensions are related with the fourier transform and not special
%   to the Zernike coefficients.  So one needs to change the 'size' of the
%   matrix that contains the pupil function as the wavelength changes (but
%   keeping the pupil size constant, so the sampling geometry is variable),
%   so that the scale on the psf is the same for each wavelength.
% 
%   What I've done by hardcoding the 550 is to say
%   that I will adjust the scale of all the psfs so that they match the
%   scale you would get at 550nm.  So yes, the 550nm is not special.  It
%   doesn't really matter whether this would be 550nm or something else
%   (or made into 'reflamda')- the only result is that  you would end up
%   with different overall scale if you used a different wavelength here.
%   So to make everything always come out with the same scale (pix size on
%   the psf) it seemed good to me to just pick a wavelength here to use as
%   the default for rescaling when computing the other wavelengths. This
%   part only affects the scale of the output and not any other part of
%   the computation. [Note from DB, I changed the hard coded 550 to setScaleWl,
%   and set that to 550.]
% 
%   At some point I tried making so that the user could put the scale
%   (arcminperpix) and the number of points across the pupil as an input
%   parameters and have the code work out what it need to do to achieve
%   this, but I ended up having a hard time finding a way to do this that
%   worked right all the time (the psfs weren't aligned all the time, I
%   guess I way having rounding issues or something).
setScaleWl = 550;
arcminperpix = (180*60/3.1416)*setScaleWl*.001*.001/sizeOfFieldMM; 
psf = zeros(sizeOfFieldPixels,sizeOfFieldPixels,length(wls));
areapix = zeros(length(wls));
areapixapod = zeros(length(wls));
for wl = 1:length(wls)
    % Get the pupil function
    zcoeffs(4) = doriginal + microns(wl);                      % Add in the appropriate LCA to the initial zernike defocus
	tempSizeOfFieldMM = sizeOfFieldMM*wls(wl)/setScaleWl;      % Rescaling so that PSF pixel dimension is constant with wavelength.
    [pupilfunc,areapix(wl),areapixapod(wl)] = ComputePupilFunctionFromZernike(zcoeffs,measpupilMM,calcpupilMM,wls(wl),sizeOfFieldPixels,tempSizeOfFieldMM,sceParams);

    % Convert to psf
    amp=fft2(pupilfunc);
    int=(amp .* conj(amp));
    psf(:,:,wl) = real(fftshift(int));
    
    % Get strehl ratio, based on Heidi's code
    strehl(wl) = max(max(psf(:,:,wl)./(areapixapod^2)));
    
    % Scale so that psf sums to unity before SCE.  
    % The variable sceFrac tells you how much light
    % is lost if you correct for the SCE.
    sceFrac(wl) = (areapixapod(wl)/areapix(wl));
    psf(:,:,wl) = psf(:,:,wl)/sum(sum(psf(:,:,wl)));
end

end


