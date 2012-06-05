function wvf = wvfComputePSF(wvf)
% wvr = wvfComputePSF(wvf)
%
% Compute the psf for the wvf object.   If the psf is already
% computed and not stale, this will return fast.  Otherwise it computes and
% stores.
%
% The point spread function is computed for each of the wavelengths listed
% in the input wvf structure. The PSF computation is based on 10 orders of
% Zernike coefficients specified to the OSA standard.
%
% See also wvfGet, wvfCreate, wvfSet, wvfComputePupilFunction
%
% Based on code provided by Heidi Hofer.
%
% 8/20/11 dhb      Rename function and pull out of supplied routine. Reformat comments.
% 9/5/11  dhb      Rename.  Rewrite for wvf i/o.
% 6/2/12  dhb      Simplify greatly given new get/set conventions.
%
% Copyright Wavefront Toolbox Team 2012

% Only do this if we need to -- it might already be computed and stored
if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE || ~isfield(wvf,'pupilfunc') || ~isfield(wvf,'PUPILFUNCTION_STALE') || wvf.PUPILFUNCTION_STALE) 
   
    % Make sure pupil function is computed.  
    wvf = wvfComputePupilFunction(wvf);
    
    wave = wvfGet(wvf,'wave');
    nWave = wvfGet(wvf,'nwave');
    psf = cell(nWave,1);
    for wl = 1:nWave
        % Convert the pupil function to the PSF  Just this simple.
        % Scale so that psf sums to unity.
        pupilfunc{wl} = wvfGet(wvf,'pupil function',wl);
        amp = fft2(pupilfunc{wl});
        inten = (amp .* conj(amp));   %intensity
        psf{wl} = real(fftshift(inten));
        psf{wl} = psf{wl}/sum(sum(psf{wl}));
    end
    
    wvf.psf = psf;
    wvf.PSF_STALE = false;
end

end


