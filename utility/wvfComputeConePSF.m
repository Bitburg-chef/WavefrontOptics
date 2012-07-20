function [conePsf,coneSceFraction] = wvfComputeConePSF(wvf)% [conePsf,coneSceFraction] = wvfComputeConePSF(wvf)%% This routine finds PSFs seen by each cone class under passed spectrum weightingSpectrum.  It gets% these by taking a weighted sum of the monochromatic PSFs, where the weights are given by% the product of the LMS spectral sensitivities and the weighting spectrum.%% The returned psfs are normalized so that they sum to unity.  If you want to figure% out the relative amount of light seen by each cone class, you need to take the% spectral sensitivities into account and also the SCE if you're using that.  The% routine does return the weighted average of the monochromatic sceFraction entries% for each cone type, to make this easier.%% Note 1 on usage.  If you actually know the hyperspectral image input, you probably don't want to use this routine.% Rather, compute the monochromatic PSFs at each wavelength and do your optical blurring in the wavelength domain,% before computing cone absorbtions.  Doing so is a more accurate model of the physics.%% You would use this routine under two circumstances.  First, you might know that the image consists only of% intensity modulations of a single relative spectrum.  In this case, you could use that spectrum here and speed% things up, since you'd only have to convolve three times (one for each cone class rather than once for each% wavelength).  This case corresponds, for example, to psychophysics where achromatic contrast is manipulated.%% Second, you might only know the unblurred LMS images and not have spectral data.  Then, this routine is useful for% providing an approximation to the blurring that will occur for each cone class.  For example, your data might% originate with a high-resolution RGB camera image, which was then used to estimate LMS values at each location.% Keep in mind that what you get in that case is only an approximation, since the actual blur depends on the% full spectral image.%% Note 2 on usage.  If you want to compute a strehl ratio quantity for the LMS psfs, the most straightforward% way is to call this routine a second time using a zcoeffs vector of all zeros.  This leads to computation% of diffraction limited monochromatic psfs that are then summed just like the specified ones.  Taking the ratios% of the peaks then gives you a fairly meaningful figure of merit.%% Note: Need to think through normalization and meaning of coneSceFraction carefully.%% See wvfComputePSF.%% 7/13/07   dhb    Made into a callable function, based on code provided by Heidi Hofer.%           dhb    Remove globals, fix case of fft, get rid of some vars we don't care about%           dhb    Don't write files here, optional plot supression.% 7/14/07   dhb    Change name a little.% 12/22/09  dhb    Return monochromatic PSFs as a cell array% 8/21/11   dhb    Update% 9/7/11    dhb    Rename.  Use wvf for i/o.% 7/20/12   dhb    Got this to run again in its modern form.%% (c) Wavefront Toolbox Team 2011, 2012%% Get wavelengthswls = wvfGet(wvf,'calc wavelengths');nWls = length(wls);%% Get weighted cone fundamentals, and% normalize each so it sums to one.conePsfInfo = wvfGet(wvf,'calc cone psf info');nCones = size(conePsfInfo.T,1);coneWeight = zeros(nCones,nWls);T = SplineCmf(conePsfInfo.S,conePsfInfo.T,wls);spdWeighting = SplineSpd(conePsfInfo.S,conePsfInfo.spdWeighting,wls);for j = 1:nCones    coneWeight(j,:) = T(j,:).*spdWeighting';end%% Get psfs for each wavelength.  This comes back as a cell% array unless there is only one wavelength. psf = wvfGet(wvf,'psf');%% Get fraction of light at each wavelength lost to scesceFraction = wvfGet(wvf,'sce fraction',wls);%% Weight up cone psfs%% Need to handle case of one wavelength separately because% this doesn't come back as a cell array.if (nWls == 1)    [m,n] = size(psf);    conePsf = zeros(m,n,nCones);    for j = 1:nCones        conePsf(:,:,j) = sceFraction*coneWeight(j)*psf;    endelse    [m,n] = size(psf{1});    conePsf = zeros(m,n,nCones);    for j = 1:nCones        for wl = 1:nWls            [m1,n1] = size(psf{wl});            if (m1 ~= m || n1 ~= n)                error('Pixel size of individual wavelength PSFs does not match');            end            conePsf(:,:,j) = conePsf(:,:,j) + sceFraction(wl)*coneWeight(j,wl)*psf{wl};        end    endend% Normalize each PSF to unit volume.for j = 1:nCones    conePsf(:,:,j) = conePsf(:,:,j)/sum(sum(conePsf(:,:,j)));end% Get sceFraction for each cone typefor j = 1:nCones    coneSceFraction(j,:) = coneWeight(j,:) .* sceFraction';end