function [conepsf,arcminperpix,defocusDiopters] = ...
    ComputeOptimizedConePSF(coneWeights,criterionFraction,wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,sizeOfFieldPixels,sizeOfFieldMM,sceParams)
% [conepsf,arcminperpix,defocusDiopters] = ...
%    ComputeOptimizedConePSF(criterion,wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,sizeOfFieldPixels,sizeOfFieldMM,sceParams)
%
% Optimize the PSF seen by the cones, given the cone sensitivities, a weighting spectral power distribution, and a criterion.  Optimization is performed on
% the defocus parameter, relative to a specified nominal focus wavelength.
%
% The coneWeights is s vector with relative weights to be applied to a figure of merit for the PSF seen by each cone class.   The total figure of merit
% is the weighted sum of the figure for each cone class alone.
%
% The figure of merit is the radius of a centered and circularly averaged version of the psf that contains the specified criterionFraction
% of the total mass of the PSF seen by each cone.  The smaller this radius, the better.
%
% 8/26/11  dhb  Wrote it.
% 8/29/11  dhb  Don't need to center or circularly average here.
%          dhb  Print warning if optimal value is at search bound.
             
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
[f,conepsf,arcminperpix] = InlineMinFunction(defocusDiopters);
if (abs(defocusDiopters) >= diopterBound)
    fprintf('WARNING: defocus found in ComputeOptimizedConePSF is at search limit of %0.1f diopters\n',diopterBound)
end

    function [f,conepsf,arcminperpix,critRadius] = InlineMinFunction(x)
        defocusDiopters = x;
        [conepsf, arcminperpix] = ...
            ComputeConePSFFromZernike(wls,T_cones,weightingSpectrum,zcoeffs,measpupilMM,calcpupilMM,nominalFocusWl,defocusDiopters,...
            sizeOfFieldPixels,sizeOfFieldMM,sceParams);
        nCones = size(T_cones,1);
        f = 0;
        for j = 1:nCones
            %temppsf = CircularlyAveragePSF(CenterPSF(conepsf(:,:,j)));
            critRadius(j) = FindPSFCriterionRadius(conepsf(:,:,j),criterionFraction);
            f = f + coneWeights(j)*critRadius(j);
        end
            
    end
end