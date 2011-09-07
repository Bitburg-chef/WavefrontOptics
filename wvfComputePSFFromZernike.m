function wvfParams = wvfComputePSF(wvfParams)
% wvfParams = wvfComputePSF(wvfParams)
%
% Computes the monochromatic psf over the calculated pupil size for 10 orders of Zernike
% coeffcients specified to the OSA standard. Includes SCE (Stiles-Crawford Effect) if specified.
%
% The psfs at each wavelength are scaled so that each sums to unity.  If SCE is specified,
% you can multiply by sceFrac at that wavelength to incorporate the fraction of light absorbed.
%
% Required input fields for wvfParams struct - see comment in wvfComputePupilFunction for more details.
%   zcoeffs -           Zernike coefficients.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%   caclpupilsize -     Size over which returned pupil function is calculated, in MM.
%   wls -               Column vector of wavelengths over which to compute, in NANOMETERS.
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   sizeOfFieldPixels - Linear size of square image over which the pupil function is computed.
%   sizeOfFieldMM -     Size of square image over which the pupile function is computed in MM.
%
% Optional input fields for wvfParams struct
%   sceParams -         Parameter structure for Stiles-Crawford correction.  If missing or set to empty,
%                       no correction and is set to empty on return.
%
% Output fields set in wvfParams struct
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
% Note: These calculations only account for lognitudinal chromatic aberration (LCA), and do not
% incoporate transverse chromatic aberration (TCA).
%
% Note: If this function is called with a zcoeffs vector that is composed of only zeros 
% the output monochromatic and polychromatic PSFs will be
% limited only by diffraction and longitudinal chromatic aberration (and the SCE if specified).
%
% Note: For better or worse, there are three ways that defocus can be
% specified to this routine.  First, the value of the 4th zernike
% coefficient is included.  Then, a value is computed based on the
% differencebetween nominalFocusWl and the wavelength for which
% the PSF is being computed. Finally, the value of defocusDiopters is 
% added in directly.  Although these are redundant, it is often most
% convenient to think in terms of one of the three ways.
%
% See also: wvfComputeConePSF, wvfComputePupilFunction, GetStilesCrawfordParamsParams, wvfGetDefocusFromWavelengthDifference
%
% Based on code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.
% 9/5/11  dhb      Rename.  Rewrite for wvfParams i/o.

%% Check optional fields
if (~isfield(wvfParams,'sceParams') || isempty(wvfParams.sceParams))
    wvfParams.sceParams = [];
end

%% Handle defocus relative to reference wavelength.
[wvfParams] = wvfGetDefocusFromWavelengthDifference(wvfParams);

% Store the original value of defocus (wvfParams.zcoeffs(4)) since it will get changed
% later when changing wavelength
doriginal = wvfParams.zcoeffs(4);

%% Get the pupil function, correcting defocus at each wavelength
%
% The value of arcminperpix is set depending on a wavelength chosen
% to set the scale.  Here is what Heidi says about this process,
% and about the renormalization of scale for each call to
% wvfComputePupilFunction:
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
arcminperpix = (180*60/3.1416)*setScaleWl*.001*.001/wvfParams.sizeOfFieldMM; 
psf = zeros(wvfParams.sizeOfFieldPixels,wvfParams.sizeOfFieldPixels,length(wvfParams.wls));
areapix = zeros(length(wvfParams.wls));
areapixapod = zeros(length(wvfParams.wls));
for wl = 1:length(wvfParams.wls)
    % Get the pupil function
    tmpWvfParams = wvfParams;
    tmpWvfParams.wls = wvfParams.wls(wl);
    tmpWvfParams.zcoeffs(4) = doriginal + wvfParams.defocusMicrons(wl);              % Add in the appropriate LCA to the initial zernike defocus
    tmpWvfParams.sizeOfFieldMM = wvfParams.sizeOfFieldMM*tmpWvfParams.wls/setScaleWl; % Rescaling so that PSF pixel dimension is constant with wavelength.
    tmpWvfParams = wvfComputePupilFunction(tmpWvfParams);
    pupilfunc(:,:,wl) = tmpWvfParams.pupilfunc;
    areapix(wl) = tmpWvfParams.areapix;
    areapixapod(wl) = tmpWvfParams.areapixapod;

    % Convert to psf
    amp=fft2(pupilfunc(:,:,wl));
    int=(amp .* conj(amp));
    psf(:,:,wl) = real(fftshift(int));
    
    % Get strehl ratio, based on Heidi's code
    strehl(wl) = max(max(psf(:,:,wl)))./(areapixapod(wl)^2);
    
    % Scale so that psf sums to unity before SCE.  
    % The variable sceFrac tells you how much light
    % is lost if you correct for the SCE.
    sceFrac(wl) = (areapixapod(wl)/areapix(wl));
    psf(:,:,wl) = psf(:,:,wl)/sum(sum(psf(:,:,wl)));
end
% setScaleWl = 550;
% arcminperpix = (180*60/3.1416)*setScaleWl*.001*.001/sizeOfFieldMM; 
% psf = zeros(sizeOfFieldPixels,sizeOfFieldPixels,length(wls));
% areapix = zeros(length(wls));
% areapixapod = zeros(length(wls));
% for wl = 1:length(wls)
%     % Get the pupil function
%     zcoeffs(4) = doriginal + microns(wl);                      % Add in the appropriate LCA to the initial zernike defocus
% 	tempSizeOfFieldMM = sizeOfFieldMM*wls(wl)/setScaleWl;      % Rescaling so that PSF pixel dimension is constant with wavelength.
%     [pupilfunc,areapix(wl),areapixapod(wl)] = ComputePupilFunction(zcoeffs,measpupilMM,calcpupilMM,wls(wl),sizeOfFieldPixels,tempSizeOfFieldMM,sceParams);
% 
%     % Convert to psf
%     amp=fft2(pupilfunc);
%     int=(amp .* conj(amp));
%     psf(:,:,wl) = real(fftshift(int));
%     
%     % Get strehl ratio, based on Heidi's code
%     strehl(wl) = max(max(psf(:,:,wl)))./(areapixapod(wl)^2);
%     
%     % Scale so that psf sums to unity before SCE.  
%     % The variable sceFrac tells you how much light
%     % is lost if you correct for the SCE.
%     sceFrac(wl) = (areapixapod(wl)/areapix(wl));
%     psf(:,:,wl) = psf(:,:,wl)/sum(sum(psf(:,:,wl)));
% end

%% Set output fields
wvfParams.psf = psf;
wvfParams.pupilfunc = pupilfunc;
wvfParams.arcminperpix = arcminperpix;
wvfParams.strehl = strehl;
wvfParams.sceFrac = sceFrac;
wvfParams.areapix = areapix;
wvfParams.areapixapod = areapixapod;

end


