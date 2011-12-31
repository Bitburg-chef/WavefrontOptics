function wvfParams = wvfGetDefocusFromWavelengthDifference(wvfParams)
% Defocus in microns, from the wavelength of nominal focus and the desired
% wavelengths and the fixed added defocusDiopters.
%
%    wvfParams = wvfGetDefocusFromWavelengthDifference(wvfParams)
%
% The two parameters nominalFocusWl and defocusDiopters are redundant, in
% the the sense that their effects get added together and you can
% accomplish the same thing with either variable.  But, sometimes it is
% more convenient to think in one way, and sometimes in the other.
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
% Examples:
%    wvfP = wvfCreate;
%    wvfP.wls = 400:10:700;
%    wvfP = wvfGetDefocusFromWavelengthDifference(wvfP)
%
% Note from DHB.  This magic code provided by Heidi Hofer.
% It does have the feature that the adjustment is 0 when
% the wavelength being evaluated matches the passed nominal
% focus wavelength.  Heidi assures me that the constants
% are correct for any pair of wavelengths.
%
% (c) Wavefront Toolbox Team 2011
%
% 8/21/11  dhb  Pulled out from code supplied by Heidi Hofer.
% 9/5/11   dhb  Rename.  Rewrite for wvfPrams i/o.

wave = wvfParams.wls;
diopters = zeros (size(wave));

% What is the reference for the value of this constant?  
constant = 1.8859 - (0.63346/(0.001*wvfParams.nominalFocusWl-0.2141));

% This appears to be a formula that calculates a defocus in diopters for
% each sample wavelength.
for ww = 1:length(wvfParams.wls)
   diopters(ww) = 1.8859 - constant - (0.63346/(0.001*wave(ww)-0.2141));  
end
% vcNewGraphWin; plot(wave,diopters); grid on

% Add in extra defocus
diopters = diopters + wvfParams.defocusDiopters;

% Convert defocus in diopters to defocus of the wavefront in microns.  We
% need a reference here to explain this.
wvfParams.defocusMicrons = zeros(size(wvfParams.wls));
wvfParams.defocusMicrons = diopters * (wvfParams.measpupilMM)^2/(16*sqrt(3));

return
