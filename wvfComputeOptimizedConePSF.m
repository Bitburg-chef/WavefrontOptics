function wvfParams = wvfComputeOptimizedConePSF(wvfParams)
% wvfParams = wvfComputeOptimizedConePSF(wvfParams)
%
% Optimize the PSF seen by the cones, given the cone sensitivities, a weighting spectral power distribution, and a criterion.  Optimization is performed on
% the defocus parameter, relative to a specified nominal focus wavelength.

% Required input fields for wvfParams struct - see comment in wvfComputePupilFunction for more details.
%   coneWeights -       vector with relative weights to be applied to a figure of merit for the
%                       PSF seen by each cone class.  The total figure of merit is the weighted sum of the
%                       figure for each cone class alone.
%   criterionFraction - The figure of merit is the radius of a centered and circularly averaged version of
%                       the psf that contains the specified criterionFraction of the total mass of the PSF seen by each cone.
%                       The smaller this radius, the better.
%   wls -               Column vector of wavelengths over which to cone sensitivity and spectrum are specified.
%   T_cones -           Cone spectral sensitivities in Psychtoolbox data format.
%   weightingSpectrum - Weighting spectrum as a column vector.
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
%   conepsf -           Calcuated psf for each cone in T_cones, third dimension indexes cone type.
%   defocusDiopters -   The defocus added in to optimize.
%   coneSceFrac -       Vector with calculated SCE fraction for each cone type.
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
% 8/26/11  dhb  Wrote it.
% 8/29/11  dhb  Don't need to center or circularly average here.
%          dhb  Print warning if optimal value is at search bound.
% 9/7/11   dhb  Rename.  Use wvfParams for i/o.

options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
if (IsCluster && matlabpool('size') > 1)
    options = optimset(options,'UseParallel','always');
end

% Initial defocus and bounds (diopters)
diopterBound = 4;
x0 = 0;
vlb = -diopterBound;
vub = -vlb;

% Optimize focus
x = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);

% Set up return values
defocusDiopters = x;
if (abs(defocusDiopters) >= diopterBound)
    fprintf('WARNING: defocus found in wvfComputeOptimizedConePSF is at search limit of %0.1f diopters\n',diopterBound)
end
[f,tmpWvfParams] = InlineMinFunction(defocusDiopters);
wvfParams = tmpWvfParams;
wvfParams.defocusDiopters = defocusDiopters;

    function [f,tmpWvfParams] = InlineMinFunction(x)
        tmpWvfParams = wvfParams;
        tmpWvfParams.defocusDiopters = x;
        tmpWvfParams = wvfComputeConePSF(tmpWvfParams);
        nCones = size(tmpWvfParams.T_cones,1);
        f = 0;
        for j = 1:nCones
            %temppsf = CircularlyAveragePSF(CenterPSF(conepsf(:,:,j)));
            critRadius(j) = FindPSFCriterionRadius(tmpWvfParams.conepsf(:,:,j),tmpWvfParams.criterionFraction);
            f = f + tmpWvfParams.coneWeights(j)*critRadius(j);
        end
        
    end
end