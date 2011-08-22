% ZernikeToolbox
%
% To do list:
%   This is in pretty clean shape, but it would be good to do a few things someday to
%   have a really complete set of optics routines.
%
%     i) Write routines for getting Williams et al PSFs and compare.
%     ii) Compare the diffraction limited PSFs from the Zernike computations to 
%     directly calculated version using PTB code.
%     iii) Implement and test an LSF to PSF routine.
%     iv) Implement the Marimont and Wandell chromatic aberration computations and compare.
%     v) Could bootstrap the precision of the average PSF estimates.
%     vi) Could try to find a functional form to produce good approximations to the PSF estimates
%     for different accommodative states, so as to have a compact and convenient way to do these
%     types of calculations in the future.
%
% 12/15/09  dhb  Wrote it as we are about to ship off the paper.
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
%   ComputePSFFromZernike -                Compute monochomatic PSFs given Zernike coeffs and other params.
%   ComputePupilFunctionFromZernike -      Compute monochromatic pupil functions given Zernike coeffs and other params.
%   GetDefocusFromWavelengthDifference -   Just as the name suggests.
%   GetStilesCrawfordParams -              Get the parameters required for incorporating Stiles-Crawford effect
%   ComputeLMSPsfs -                       Compute LMS psfs given zernike, pupil, accom wl, and weighting spectrum
%   FindMatPeak -                          Find the coordinates of the maximum in a matrix.
%   FindPSFCriterionRadius -               Find radius that includes passed fraction of PSF mass
%
% Data
%   sampleZernikeCoeffs.txt -              Zernike coefficient data for 9 subjects in OSA format, measured by Heidi Hofer.
%                                          Data for each subject, 65 coefficients in a column.  The pupil size for
%                                          the measurements was 6 mm for all subjects.  Tip, tilt, and defocus
%                                          coefficients are set to zero.
%   xZernikePsfs -                         Computed by GenerateRawPSFs.  Do it on a cluster if possible.