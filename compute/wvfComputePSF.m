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
%   sizeOfFieldPixels - Linear size of square image over which the pupil function is computed.
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
% aberration (LCA), and do not incoporate transverse chromatic aberration
% (TCA).
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
% These are redundant, it is often most convenient to think in terms of one
% of the three ways. 
%
% See also: wvfComputeConePSF, wvfComputePupilFunction, sceGetParamsParams,
% wvfGetDefocusFromWavelengthDifference 
%
% Based on code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.
% 9/5/11  dhb      Rename.  Rewrite for wvfP i/o.

%% Check optional fields
if (~isfield(wvfP,'sceParams') || isempty(wvfP.sceParams))
    wvfP.sceParams = [];
end

%% Handle defocus relative to reference wavelength.
defocusMicrons = wvfGet(wvfP,'defocus distance','um');
% wvfP.defocusMicrons = defocusMicrons; % Old
wvfP = wvfSet(wvfP,'defocus microns',defocusMicrons);

% wvfGetDefocusFromWavelengthDifference(wvfP);

% Get the original  defocus coefficient (wvfP.zcoeffs(4)).
% It will get changed later when changing wavelength
doriginal = wvfGet(wvfP,'zcoeffs',4);
% doriginal = wvfP.zcoeffs(4);

%% Get the pupil function, correcting defocus at each wavelength
%
% The value of arcminperpix is set depending on a wavelength chosen
% to set the scale.  Here is what Heidi says about this process,
% and about the renormalization of scale for each call to
% wvfComputePupilFunction:
%
%   The setting of the scale and then the change in the loop with
%   wavelength is to make sure that the scale on the psfs stays the same as
%   the wavelength is changed - since the scale is related to lamda/(linear
%   dimension of matrix) - this is just related to the way the dimensions
%   are related with the fourier transform and not special to the Zernike
%   coefficients.  So one needs to change the 'size' of the matrix that
%   contains the pupil function as the wavelength changes (but keeping the
%   pupil size constant, so the sampling geometry is variable), so that the
%   scale on the psf is the same for each wavelength.
% 
%   What I've done by hardcoding the 550 is to say that I will adjust the
%   scale of all the psfs so that they match the scale you would get at
%   550nm.  So yes, the 550nm is not special.  It doesn't really matter
%   whether this would be 550nm or something else (or made into
%   'reflamda')- the only result is that  you would end up with different
%   overall scale if you used a different wavelength here. So to make
%   everything always come out with the same scale (pix size on the psf) it
%   seemed good to me to just pick a wavelength here to use as the default
%   for rescaling when computing the other wavelengths. This part only
%   affects the scale of the output and not any other part of the
%   computation. [Note from DB, I changed the hard coded 550 to setScaleWl,
%   and set that to 550.]
% 
%   At some point I tried making so that the user could put the scale
%   (arcminperpix) and the number of points across the pupil as an input
%   parameters and have the code work out what it need to do to achieve
%   this, but I ended up having a hard time finding a way to do this that
%   worked right all the time (the psfs weren't aligned all the time, I
%   guess I was having rounding issues or something).

setScaleWl = 550;   % Nanometers

% It looks like this has three parts
%   (180/pi)*60   setScaleWl -> millimeters      1/fieldMM
%
% Here is my understanding
%   (180/pi) is 1 deg in radians, 60 gets us minutes in radians.
%   (deg/rad)*min/deg ->  min/rad
%
%   setScaleWl*10^-6 is nm * mm/nm = mm
%   (1/mm) is the total field of view
%
%   So we have min/rad units when we are done.  I guess if 1 minute is 1
%   pixel, then we have arc minutes per pixel.
%
% This is probably close, but not quite right.  Maybe we can fix this
% comment.
% arcminperpix = (180*60/3.1416)*setScaleWl*.001*.001/wvfP.sizeOfFieldMM; 

wave    = wvfGet(wvfP,'wave');
nWave   = wvfGet(wvfP,'nwave');
nPixels = wvfGet(wvfP,'npixels');

psf     = zeros(nPixels,nPixels,nWave);
areapix = zeros(nWave,1);
areapixapod = zeros(nWave,1);
strehl  = zeros(nWave,1);
sceFrac = zeros(nWave,1);

for wl = 1:nWave
    
    % Get the pupil function.  
    
    % Set up a tmp parameter structure so we can do it one wave at a time.
    tmpWvfParams = wvfP;
    tmpWvfParams = wvfSet(tmpWvfParams,'wave',wave(wl));
    
    % Add in the appropriate LCA to the initial zernike defocus
    % Don't understand the following.
    tmpWvfParams.zcoeffs(4) = doriginal + wvfP.defocusMicrons(wl); 
    
    % Rescaling so that PSF pixel dimension is constant with wavelength.
    tmpWvfParams.sizeOfFieldMM = wvfP.sizeOfFieldMM*tmpWvfParams.wls/setScaleWl;
    
    % Compute the pupil function
    tmpWvfParams = wvfComputePupilFunction(tmpWvfParams);
    pupilfunc(:,:,wl) = wvfGet(tmpWvfParams,'pupil function');
    areapix(wl)       = wvfGet(tmpWvfParams,'area pix');
    areapixapod(wl)   = wvfGet(tmpWvfParams,'area pixapod');

    % Convert to the pupil function to the PSF  Just this simple.
    amp = fft2(pupilfunc(:,:,wl));
    inten = (amp .* conj(amp));   %intensity
    psf(:,:,wl) = real(fftshift(inten));
    
    % Get strehl ratio, based on Heidi's code
    strehl(wl) = max(max(psf(:,:,wl)))./(areapixapod(wl)^2);
    
    % Scale so that psf sums to unity before SCE.  
    % The variable sceFrac tells you how much light
    % is lost if you correct for the SCE.
    sceFrac(wl) = (areapixapod(wl)/areapix(wl));
    psf(:,:,wl) = psf(:,:,wl)/sum(sum(psf(:,:,wl)));
end

%% Set output fields
wvfP = wvfSet(wvfP,'pupil function',pupilfunc);
wvfP = wvfSet(wvfP,'psf',psf);
wvfP = wvfSet(wvfP,'strehl',strehl);

% Derived, not needed? wvfP.arcminperpix = arcminperpix;
% wvfP.strehl = strehl;
% wvfP.sceFrac = sceFrac;
% Don't understand yet. (BW)
% wvfP.areapix = areapix;
% wvfP.areapixapod = areapixapod;

end


