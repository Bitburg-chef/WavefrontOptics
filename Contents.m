% WavefrontOpticsToolbox
%
% Programs (demonstrate use)
%   s_wvfComputeAverageObserverConePSF -      Compute optimized cone psf for an average observer. 
%   s_wvfComputeConePSFFromZernikeTest -      Test the cone/spectrum weighted psf routines.
%   s_wvfComputePSFFromZernikeTest -          Basic test that underlying Zernike monochromatic psf routines behave sensibly.
%
% Wavefront optics suite
%   wvfComputeConePSFFromZernike -            Compute cone psfs given zernike, pupil, accom wl, and weighting spectrum
%   wvfComputeOptimizedConePSF -              Compute focus optimized cone psfs given zernike, pupil, accom wl, and weighting spectrum
%   wvfComputePSFFromZernike -                Compute monochomatic PSFs given Zernike coeffs and other params.
%   wvfComputePupilFunctionFromZernike -      Compute monochromatic pupil functions given Zernike coeffs and other params.
%   wvfGetDefocusFromWavelengthDifference -   Just as the name suggests.
%
% Helper routines (in subdir helpers)
%   AverageOpticalPSFs -                      Average passed PSFs
%   CenterPSF -                               Put max of PSF at center of image
%   CircularlyAveragePSF -                    Compute circular average of a PSF
%   FindMatPeak -                             Find the coordinates of the maximum in a matrix.
%   FindPSFCriterionRadius -                  Find radius that includes passed fraction of PSF mass
%   GetStilesCrawfordParamsParams -           Get the parameters required for incorporating Stiles-Crawford effect
%
% Data
%   sampleZernikeCoeffs.txt -                 Zernike coefficient data for 9 subjects in OSA format, measured by Heidi Hofer.
%                                             Data for each subject, 65 coefficients in a column.  The pupil size for
%                                             the measurements was 6 mm for all subjects.  Tip, tilt, and defocus
%                                             coefficients are set to zero.
