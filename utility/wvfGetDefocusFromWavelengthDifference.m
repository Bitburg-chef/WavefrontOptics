function defocusMicrons = wvfGetDefocusFromWavelengthDifference(wvfP)
% Defocus in microns, from the wavelength of nominal focus and the desired
% wavelengths and the fixed added defocusDiopters.
%
%    defocusMicrons = wvfGetDefocusFromWavelengthDifference(wvfP)
%
% The two parameters nominalFocusWl and defocusDiopters are redundant, in
% the the sense that their effects get added together and you can
% accomplish the same thing with either variable.  But, sometimes it is
% more convenient to think in one way, and sometimes in the other.
%
% The conversion to microns is because these are the units we assume that
% the pupil function is measured in.
%
% Required input fiels for wvfP struct - see comment in wvfComputePupilFunction for more details.
%   wls -               Column vector of wavelengths (in nm) to compute
%   nominalFocusWl -    Wavelength (in nm) of nominal focus.
%   defocusDiopters -   Defocus to add in (signed), in diopters.
%   measpupilMM -       Size of pupil characterized by the coefficients, in MM.
%
% Output fields set in wvfP struct
%   defocusMicrons -    Defocus added in to zcoeffs(4) at each wavelength, in microns.
%
% Examples:
%    wvfP = wvfCreate;
%    wvfP = wvfSet(wvfP,'wave',400:10:700);
%    wvfP = wvfGetDefocusFromWavelengthDifference(wvfP)
%
% Note from DHB.  This magic code provided by Heidi Hofer. It does have the
% feature that the adjustment is 0 when the wavelength being evaluated
% matches the passed nominal focus wavelength.  Heidi assures me that the
% constants are correct for any pair of wavelengths.
%
% 8/21/11  dhb  Pulled out from code supplied by Heidi Hofer.
% 9/5/11   dhb  Rename.  Rewrite for wvfPrams i/o.
%
% (c) Wavefront Toolbox Team 2011

wave = wvfGet(wvfP,'wave'); 
diopters = zeros (size(wave));

% What is the reference for the value of this constant?  
nominalFocusWl = wvfGet(wvfP,'infocuswave');
constant = 1.8859 - (0.63346/(0.001*nominalFocusWl-0.2141));

% This appears to be a formula that calculates a defocus in diopters for
% each sample wavelength.
for ww = 1:length(wave)
   diopters(ww) = 1.8859 - constant - (0.63346/(0.001*wave(ww)-0.2141));  
end
% vcNewGraphWin; plot(wave,diopters); grid on

% Add in extra defocus from using non-nominal focus wavelength 
% to the original wvfP-specified defocus 
defocusDiopters = wvfGet(wvfP,'defocusdiopters');
diopters = diopters + defocusDiopters;

% Convert defocus in diopters to defocus of the wavefront in microns.
wvfP.defocusMicrons = zeros(size(wave)); %this looks unused (KP)

pupilMM = wvfGet(wvfP,'measured pupil','mm');
defocusMicrons = diopters * (pupilMM)^2/(16*sqrt(3));

% References for equation:
% http://www.telescope-optics.net/monochromatic_eye_aberrations.htm
% (accessed 3/11/12) 
% http://www.journalofvision.org/content/10/5/4.full (accessed 3/11/12)
% (under Zernike defocus eqn 1)
% Several other references, like Thibos, have referred to the same equation
% but I have not found a full derivation.
% http://research.opt.indiana.edu/Library/HVO/Handbook.html
% (accessed 3/11/12) comes close in the "numerical example: defocus"
% section, but is still missing a factor of 2 in the denominator.
% (KP)

return
