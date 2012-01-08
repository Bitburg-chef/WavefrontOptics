% WavefrontOpticsToolbox
%
% Programs (demonstrate use)
%   s_wvfComputeAverageObserverConePSF -      Compute optimized cone psf for an average observer. 
%   s_wvfComputeConePSFTest -                 Test the cone/spectrum weighted psf routines.
%   s_wvfComputePSFTest -                     Basic test that underlying Zernike monochromatic psf routines behave sensibly.
%
% Wavefront (wvf) optics suite
%   wvfComputeConePSF -                       Compute cone psfs given Zernike coeffs, pupil, accom wl, and weighting spectrum
%   wvfComputeOptimizedConePSF -              Compute focus optimized cone psfs given Zernike coeffs, pupil, accom wl, and weighting spectrum
%   wvfComputeOptimizedPSF -                  Compute focus optimized (at a specified wl) monochromatic PSFs, ginve Zernike coeffs, etc.
%   wvfComputePSF -                           Compute monochomatic PSFs given Zernike coeffs, etc.
%   wvfComputePupilFunction -                 Compute monochromatic pupil functions given Zernike coeffs, etc.
%   wvfGetDefocusFromWavelengthDifference -   Just as the name suggests.
%
% Stiles-Crawford Effect (sce)
%   sceGetParamsParams -                      Get the parameters required for incorporating Stiles-Crawford effect
%
% Operate on point spread functions (psf)
%   psfAverageMultiple -                      Average passed PSFs
%   psfCenter -                               Put max of PSF at center of image
%   psfCircularlyAverage -                    Compute circular average of a PSF
%   psfFindCriterionRadius -                  Find radius that includes passed fraction of PSF mass
%   psfFindPeak -                             Find the coordinates of the maximum in a matrix.
%
% Data
%   autrusseauStandardObserver.txt -          Zernike coefficient data for the "standard observer" of Autrusseau, Thibos, & Shevell (2011),
%                                             VisionResearch, 51, pp. 2282-2294, Table 1.  These seem to have defocus, tip, and tilt included
%                                             so it is worth some checking to make sure the convention matches the code here.
%   sampleZernikeCoeffs.txt -                 Zernike coefficient data for 9 subjects in OSA format, measured by Heidi Hofer.
%                                             Data for each subject, 65 coefficients in a column.  The pupil size for
%                                             the measurements was 6 mm for all subjects.  Tip, tilt, and defocus
%                                             coefficients are set to zero.
