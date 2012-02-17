function [outPSF,peakRow,peakCol] = psfCenter(inPSF)
% Put the maximum of a PSF at the center of the two D grid.
%
%   outPSF = psfCenter(inPSF)
%
% There should be an inverse to this.  The extrapolated values are set to
% 0. 
%
% 8/26/07  dhb  Wrote it.
% 8/22/11  dhb  A 'round' should be a 'floor', I think.
%
% (c) Wavefront Toolbox Team, 2012

% Use interpolation to recenter
[peakRow,peakCol] = psfFindPeak(inPSF);
% vcNewGraphWin; mesh(inPSF)

% Interpolate data so peak is at near 0,0.  Extrapolated values are assumed
% to be 0.
[m,n] = size(inPSF);
xIn = ((1:n)-peakCol);
yIn = ((1:m)-peakRow);
xOut = ((1:n)-(floor(n/2)+1));
yOut = ((1:m)-(floor(m/2)+1));

outPSF = interp2(xIn,yIn',inPSF,xOut,yOut','linear',0);

return
