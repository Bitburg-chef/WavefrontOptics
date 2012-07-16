function wvf = wvfComputeConePSF(wvf)% wvf = wvfComputeConePSF(wvf)%% This routine finds PSFs seen by each cone class under passed spectrum weightingSpectrum.  It gets% these by taking a weighted sum of the monochromatic PSFs, where the weights are given by% the product of the LMS spectral sensitivities and the weighting spectrum.%% The returned psfs are normalized so that they sum to unity.  If you want to figure% out the relative amount of light seen by each cone class, you need to take the% spectral sensitivities into account and also the SCE if you're using that.  The% routine does return the weighted average of the monochromatic sceFrac entries% for each cone type, to make this easier.%% Note 1 on usage.  If you actually know the hyperspectral image input, you probably don't want to use this routine.% Rather, compute the monochromatic PSFs at each wavelength and do your optical blurring in the wavelength domain,% before computing cone absorbtions.  Doing so is a more accurate model of the physics.%% You would use this routine under two circumstances.  First, you might know that the image consists only of% intensity modulations of a single relative spectrum.  In this case, you could use that spectrum here and speed% things up, since you'd only have to convolve three times (one for each cone class rather than once for each% wavelength).  This case corresponds, for example, to psychophysics where achromatic contrast is manipulated.%% Second, you might only know the unblurred LMS images and not have spectral data.  Then, this routine is useful for% providing an approximation to the blurring that will occur for each cone class.  For example, your data might% originate with a high-resolution RGB camera image, which was then used to estimate LMS values at each location.% Keep in mind that what you get in that case is only an approximation, since the actual blur depends on the% full spectral image.%% Note 2 on usage.  If you want to compute a strehl ratio quantity for the LMS psfs, the most straightforward% way is to call this routine a second time using a zcoeffs vector of all zeros.  This leads to computation% of diffraction limited monochromatic psfs that are then summed just like the specified ones.  Taking the ratios% of the peaks then gives you a fairly meaningful figure of merit.%% See wvfComputePSF.%% 7/13/07   dhb    Made into a callable function, based on code provided by Heidi Hofer.%           dhb    Remove globals, fix case of fft, get rid of some vars we don't care about%           dhb    Don't write files here, optional plot supression.% 7/14/07   dhb    Change name a little.% 12/22/09  dhb    Return monochromatic PSFs as a cell array% 8/21/11   dhb    Update% 9/7/11    dhb    Rename.  Use wvf for i/o.%% (c) Wavefront Toolbox Team 2011, 2012% Get wavelengthswls = wvfGet(wvf,'calc wavelengths');nWls = length(wls);% Get weighted cone fundamentals, and% normalize each so it sums to one.nCones = size(wvf.T_cones,1);coneWeight = zeros(nCones,nWls);for j = 1:nCones    coneWeight(j,:) = wvf.T_cones(j,:).*wvf.weightingSpectrum';endfor j = 1:nCones    coneWeight(j,:) = coneWeight(j,:)./sum(coneWeight(j,:));end% Get psfs for each cone typewvf = wvfComputePSF(wvf);psf = wvfGet(wvf,'psf');[m,n] = size(psf(:,:,1));wvf.conepsf = zeros(m,n,nCones);for j = 1:nCones    for wl = 1:nWls        wvf.conepsf(:,:,j) = wvf.conepsf(:,:,j) + wvf.sceFrac(wl)*coneWeight(j,wl)*wvf.psf(:,:,wl);    endend% Normalize.  Because of the SCE, they might not be.for j = 1:nCones    wvf.conepsf(:,:,j) = wvf.conepsf(:,:,j)/sum(sum(wvf.conepsf(:,:,j)));end% Get sceFrac for each cone typefor j = 1:nCones    wvf.coneSceFrac(j,:) = coneWeight(j,:) .* wvf.sceFrac;end