function outPSF = psfCenter(inPSF)
% outPSF = psfCenter(inPSF)
%
% Put the maximum of a PSF at the center of the two D grid.
%
% 8/26/07  dhb  Wrote it.
% 8/22/11  dhb  A 'round' should be a 'floor', I think.

% Use interpolation to recenter
[peakRow,peakCol] = psfFindPeak(inPSF);
[m,n] = size(inPSF);
xIn = ((1:n)-peakCol);
yIn = ((1:m)-peakRow);
xOut = ((1:n)-(floor(n/2)+1));
yOut = ((1:m)-(floor(m/2)+1));
outPSF = interp2(xIn,yIn',inPSF,xOut,yOut','linear',0);
