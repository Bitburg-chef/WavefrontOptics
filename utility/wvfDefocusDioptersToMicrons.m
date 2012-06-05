function microns = wvfDefocusDioptersToMicrons(diopters,pupilSizeMM)
% microns = wvfDefocusDioptersToMicrons(diopters,pupilSizeMM)
%
% Convert defocus expressed in microns to defocus expressed in microns.
% The latter is suitable for adding into the 4th Zernike coefficient.
%
% The pupil size should be that used to normalize the radius of the
% Zernike coefficients, that is the size with respect to which the
% meausurements were made.
%
% 6/5/12  dhb  Wrote it as separate function.

microns = diopters*(pupilSizeMM )^2/(16*sqrt(3));
