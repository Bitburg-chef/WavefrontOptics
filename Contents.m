% WavefrontOpticsToolbox
%
% Calling routines
%   ComputeLMSPSFFromZernikeTest -         Test the LMS psf routines.
%   ComputePSFFromZernikeTest -            Basic test that underlying Zernike monochromatic psf routines behave sensibly.
%
% Helper routines
%   AverageOpticalPSFs -                   Average passed PSFs
%   CenterPSF -                            Put max of PSF at center of image
%   CircularlyAveragePSF -                 Compute circular average of a PSF
%   ComputeLMSPSFFromZernike -             Compute LMS psfs given zernike, pupil, accom wl, and weighting spectrum
%   ComputePSFFromZernike -                Compute monochomatic PSFs given Zernike coeffs and other params.
%   ComputePupilFunctionFromZernike -      Compute monochromatic pupil functions given Zernike coeffs and other params.
%   GetDefocusFromWavelengthDifference -   Just as the name suggests.
%   GetStilesCrawfordParams -              Get the parameters required for incorporating Stiles-Crawford effect
%   FindMatPeak -                          Find the coordinates of the maximum in a matrix.
%   FindPSFCriterionRadius -               Find radius that includes passed fraction of PSF mass
%
% Data
%   sampleZernikeCoeffs.txt -              Zernike coefficient data for 9 subjects in OSA format, measured by Heidi Hofer.
%                                          Data for each subject, 65 coefficients in a column.  The pupil size for
%                                          the measurements was 6 mm for all subjects.  Tip, tilt, and defocus
%                                          coefficients are set to zero.
