function wvfParams = wvfGetDefocusFromWavelengthDifference(wvfParams)
% wvfParams = wvfGetDefocusFromWavelengthDifference(wvfParams)
%
% Get the defocus in microns, from the wavelength of nominal focus
% and the desired wavelengths and the fixed added defocusDiopters.
%
% The two parameters nominalFocusWl and defocusDiopters are redundant, in the
% the sense that their effects get added together and you can accomplish
% the same thing with either variable.  But, sometimes it is more
% convenient to think in one way, and sometimes in the other.
%
% The conversion to microns is because these are the units we assume that
% the pupil function is measured in.
%
% Required input fiels for wvfParams struct - see comment in wvfComputePupilFunction for more details.
%   wls -               Column vector of wavelengths over which to compute, in NANOMETERS.
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%
% Output fields set in wvfParams struct
%   defocusMicrons -    Defocus added in to zcoeffs(4) at each wavelength, in microns.
%
% 8/21/11  dhb  Pulled out from code supplied by Heidi Hofer.
% 9/5/11   dhb  Rename.  Rewrite for wvfPrams i/o.

% Note from DHB.  This magic code provided by Heidi Hofer.
% It does have the feature that the adjustment is 0 when
% the wavelength being evaluated matches the passed nominal
% focus wavelength.  Heidi assures me that the constants
% are correct for any pair of wavelengths.
diopters = zeros (size(wvfParams.wls));
constant = 1.8859-(0.63346/(0.001*wvfParams.nominalFocusWl-0.2141));
for wl = 1:length(wvfParams.wls)
   diopters(wl) = 1.8859 - constant - (0.63346/(0.001*wvfParams.wls(wl)-0.2141));  
end

% Add in extra
diopters = diopters + wvfParams.defocusDiopters;

% Convert defocus diopters to microns
wvfParams.defocusMicrons = zeros(size(wvfParams.wls));
wvfParams.defocusMicrons = diopters * (wvfParams.measpupilMM*wvfParams.measpupilMM)/(16*sqrt(3));
