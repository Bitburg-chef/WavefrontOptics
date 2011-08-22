function [microns,diopters] = GetDefocusFromWavelenthDifference(wls,nominalFocusWl,defocusDiopters,measpupilMM)
% [microns,diopters] = GetDefocusFromWavelenthDifference(wls,nominalFocusWl,measpupilMM)
%
% Get the defocus in microns/diopters, from the wavelength of nominal focus
% and the desired wavelengths.
% 
% Also adds in an extra (signed) defocus passed in diopters.
%
% The conversion to microns is because these are the units we assume that
% the pupil function is measured in.
%
% Wavelengths passed in NM.
%
% 8/21/11  dhb  Pulled out from code supplied by Heidi Hofer.

% Note from DHB.  This magic code provided by Heidi Hofer.
% It does have the feature that the adjustment is 0 when
% the wavelength being evaluated matches the passed nominal
% focus wavelength.
diopters = zeros (size(wls));
constant = 1.8859-(0.63346/(0.001*nominalFocusWl-0.2141));
for wl = 1:length(wls)
   diopters(wl) = 1.8859 - constant - (0.63346/(0.001*wls(wl)-0.2141));  
end

% Add in extra
diopters = diopters + defocusDiopters;

% Convert defocus diopters to microns
microns = zeros(size(wls));
microns = diopters * (measpupilMM*measpupilMM)/(16*sqrt(3));
