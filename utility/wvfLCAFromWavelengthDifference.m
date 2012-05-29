function lcaDiopters = wvfLCAFromWavelengthDifference(wl1NM,wl2NM)
% lcaDiopters = wvfLCAFromWavelengthDifference(wl1NM,wl2NM)
%
% Longitudinal chromatic aberration (LCA), expressed in diopters, between
% two wavelengths.
%
% If the image is in focus at wl1NM, add the answer to bring it into focus at
% wl2NM.
%
% Either input argument may be a vector, but if both are vectors they need
% to have the same dimensions.
%
% Note from DHB.  This magic code provided by Heidi Hofer. It does have the
% feature that the adjustment is 0 when the wavelength being evaluated
% matches the passed nominal focus wavelength.  Heidi assures me that the
% constants are correct for any pair of wavelengths.
%
% 8/21/11  dhb  Pulled out from code supplied by Heidi Hofer.
% 9/5/11   dhb  Rename.  Rewrite for wvfPrams i/o.
% 5/29/12  dhb  Pulled out just the bit that computes diopters from wavelengths.
%
% (c) Wavefront Toolbox Team 2011

wave = wvfGet(wvfP,'wave'); 
diopters = zeros (size(wave));

return
