% WavefrontOpticsToolbox
%
% TUTORIAL SCRIPTS (subdir tutorial).  These help you learn about the code and the
% ideas behind it.
%  tutorial/t_wvfZernike -                    Basic wavefront optics ideas including Zernike expansion to describe wavefront aberrations
%
% VALIDATION SCRIPTS (subdir validate).  These test that various pieces of the toolboxes
% do the right thing.
%   validate/v_wvfDiffractionPSF -            Compare PSFs obtained with wavefront calcs to direct computation of diffraction limited PSFs
%   validate/v_wvfSpatialSampling -           Compare PSFs with different spatial sampling parameters.
%   validate/v_wvfSVNVer121TestData -         Compare PSFs computed by toolbox now to those computed long ago and far away.
%
% BASIC OPERATIONS (subdir wvf).  These are the core functions that know about the internals of the
% wvf structure.  User code should generally use these and not rely on the specific
% fieds within the structure.
%
% The get of pupil function and psf require that they first be explicitly computed.
% Note, however, that the ComputePSF routine will force a compute of the pupil function
% if has not already been computed.
%   wvf/wvfCreate -                           Create wavefront optics (wvf) structure
%   wvf/wvfComputePupilFunction -             Compute the pupil function for current structure values.
%   wvf/wvfComputePSF -                       Compute the psf for current structure values.
%   wvf/wvfGet -                              Get value from wvf structure
%   wvf/wvfPlot -                             Various useful plots of stuff in the wvf structure.
%   wvf/wvfPrint -                            Print values in wvf structure
%   wvf/wvfSet -                              Set value in wvf structure
%
% UTILITY ROUTINES (subdir utility).  These perform various little useful
% computations.  Some are passed the wvf structure, but none set its fields
% directly or access its fields other than via wvfGet.
%   utility/wvfComputeConePSF -               Compute the PSF seen by cones.
%   utility/wvfDefocusDioptersToMicrons -     Convert diopters to microns for addition into j=4 Zernike coefficient.
%   utility/wvfLCAFromWavelengthDifference -  Compute longitudinal chromatic aberration (LCA) in diopters, from wavelength difference.
%   utility/wvfWave2idx -                     Convert wvf wavelengths to indices.
%
% PSF MANIPULATION AND ANALYSIS (subdir psf).  Perform operations on psfs.  Not directly tied
% to the wvf structure.
%   psf/psfAverageMultiple -                  Find the average of multiple psfs.
%   psf/psfCenter -                           Put maximum of psf at matrix center.
%   psf/psfCircularlyAverage -                Circularly average a psf.
%   psf/psfFindCriterionRadius -              Find the radius of circularly symmetric psf that contains a specfied fraction of the psf max.
%   psf/psfFindPeak -                         Find the location of the maximum of a psf.
%
% STILES-CRAWFORD EFFECT (subdir stilescrawford).  As the name says.
%   sceCreate -                               Create Stiles-Crawford structure with various options as to data.
%   sceGet -                                  Get info from Stiles-Crawford structure.

% (c) Wavefront Toolbox Team 2011, 2012
%
% ****
% STUFF BELOW HERE NEEDS UPDATING.
%
% Programs (demonstrate use)
%   s_wvfComputeAverageObserverConePSF -      Compute optimized cone psf for an average observer. 
%   s_wvfComputeConePSFTest -                 Test the cone/spectrum weighted psf routines.
%   s_wvfComputePSFTest -                     Basic test that underlying Zernike monochromatic psf routines behave sensibly.
%
% Wavefront (wvf) optics suite
%   wvfComputeOptimizedConePSF -              Compute focus optimized cone psfs given Zernike coeffs, pupil, accom wl, and weighting spectrum
%   wvfComputeOptimizedPSF -                  Compute focus optimized (at a specified wl) monochromatic PSFs, ginve Zernike coeffs, etc.
%
% Data
%   autrusseauStandardObserver.txt -          Zernike coefficient data for the "standard observer" of Autrusseau, Thibos, & Shevell (2011),
%                                             VisionResearch, 51, pp. 2282-2294, Table 1.  These seem to have defocus, tip, and tilt included
%                                             so it is worth some checking to make sure the convention matches the code here.  Could their
%                                             coeficient 0 correspond to our 3?
%   sampleZernikeCoeffs.txt -                 Zernike coefficient data for 9 subjects in OSA format, measured by Heidi Hofer.
%                                             Data for each subject, 65 coefficients in a column.  The pupil size for
%                                             the measurements was 6 mm for all subjects.  Tip, tilt, and defocus
%                                             coefficients are set to zero.
