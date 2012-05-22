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

%% Handle defocus relative to reference wavelength.
defocusMicrons = wvfGet(wvfP,'defocus distance','um');
% This is weird.
wvfP = wvfSet(wvfP,'defocus microns',defocusMicrons);

% wvfGetDefocusFromWavelengthDifference(wvfP)
% is used in wvfGet(wvfP, 'defocus distance','um') to calculate
% defocusMicrons

% Get the original  defocus coefficient (wvfP.zcoeffs(4)).
% It will get changed later when changing wavelength
doriginal = wvfGet(wvfP,'zcoeffs',4);
% doriginal = wvfP.zcoeffs(4);

%% Get the pupil function, correcting defocus at each wavelength
%
% These comments no longer belong here.  I left them because this is where
% they started. But the key stuff is now moved into wvfGet in the
% 'arcminperpix' section.
%
% The value of arcminperpix is set depending on a wavelength chosen to set
% the scale.  Here is what Heidi says about this process, and about the
% renormalization of scale for each call to wvfComputePupilFunction:
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
wave    = wvfGet(wvfP,'wave');
nWave   = wvfGet(wvfP,'nwave');
% nPixels = wvfGet(wvfP,'npixels');

psf         = cell(nWave,1);
pupilfunc   = cell(nWave,1);

areapix = zeros(nWave,1);
areapixapod = zeros(nWave,1);
strehl  = zeros(nWave,1);
sceFrac = zeros(nWave,1);

dm = wvfGet(wvfP,'defocus microns');

% We will set tmpWvfParams to a single wavelength and loop.
w = waitbar(0,'Pupil functions');
for wl = 1:nWave
    
    waitbar(wl/nWave,w,sprintf('Pupil function %d',wave(wl)));
    
    % Set up a tmp structure so we can loop on wavelength
    tmpWvfParams = wvfSet(wvfP,'wave',wave(wl));
    
    % Add in the appropriate LCA to the initial zernike defocus
   
    % 3/10/12 KP: this looks like they are adjusting the defocus term in
    % the zernike coefficients (term #4) by adding Longitudinal Chromatic 
    % Aberration.
    % Takes the original Zernike coeff of "Defocus" and adds the additional 
    % defocus caused by measuring the PSF for the specified wavelength (wl)
    % rather than the given nominal focus wavelength (nominalFocuswl) using 
    % wvfGetDefocusFromWavelengthDifference (which also adds in
    % contribution from "defocus diopter" value in wvfP)
    tmpWvfParams.zcoeffs(4) = doriginal + dm(wl);
    % wvfP.defocusMicrons(wl); 
    
    % Rescaling so that PSF pixel dimension is constant with wavelength.
    % OK, this needs a much better explanation and justification.  It
    % probably breaks some of the wavelength-dependent calculations - BW.   
    % npixels = wvfGet(wvfP,'npixels');
    
    % ***** START HERE ******
    % I think the calculation should always be part of the get.  We
    % shouldn't have to adjust it on the fly like this.
    % So, 
    %    wvfGet(wvfP,'field sample size mm',waveIdx)
    % should return what is below.  Then the next line would be irrelevant.
    % newFieldSampleSize = wvfGet(wvfP,'field sample size mm')*wave(wl)/setScaleWl;
    % newFieldSampleSize = wvfGet(wvfP,'field sample size mm');
    % tmpWvfParams       = wvfSet(tmpWvfParams,'fieldSampleSize',newFieldSampleSize);
    %
    % newFieldSizeMM = wvfGet(tmpWvfParams,'field sample size mm')*npixels;
    % tmpWvfParams   = wvfSet(tmpWvfParams,'field size mm',newFieldSizeMM);
    %
    % Compute the pupil function
    tmpWvfParams      = wvfComputePupilFunction(tmpWvfParams);
    pupilfunc{wl}     = wvfGet(tmpWvfParams,'pupil function',1); 
%     areapix(wl)       = wvfGet(tmpWvfParams,'area pix');
%     areapixapod(wl)   = wvfGet(tmpWvfParams,'area pixapod');

    % Convert to the pupil function to the PSF  Just this simple.
    % Scale so that psf sums to unity before SCE is applied.
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
wvfP = wvfSet(wvfP,'pupil function',pupilfunc);
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


