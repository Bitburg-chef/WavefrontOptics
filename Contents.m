% WavefrontOpticsToolbox
%
% To do list:
%
%     i) Write routines for getting Williams et al PSFs and compare.
%     iv) Compare what these generate to Marimont and Wandell chromatic aberration computations.
%     vi) Could try to find a functional form to produce good approximations to the PSF estimates
%     for different accommodative states, so as to have a compact and convenient way to do these
%     types of calculations in the future.
%
% Calling routines
%   ComputePSFFromZernikeBasicTest -       Basic test that underlying Zernike routines behave sensibly.
%   GenerateRawPSFs -                      Generate big files of PSFs for observers/focusWls from zernike data
%   Get552MonoPSF -                        Get the average PSF for 552 nm in best focus, compare to literature
%   GetLMSPSFs -                           Get average PSFs seen by L, M, and S cones for various accommodation optimzations
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
