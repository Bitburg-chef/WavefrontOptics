function outPSF = CircularlyAveragePSF(inPSF)
% outPSF = CircularlyAveragePSF(inPSF)
%
% As the name suggests.
%
% 7/19/07   dhb  Wrote it.
% 12/22/09  dhb  Fix bug in how peakRow and peakCol are computed.
% 12/22/09  dhb  Make computation a little more fine grained.

% Define quantization.  Four was used in early code, but 1 makes more sense.
quantizationFactor = 1;

% Make a circularly symmetric version of average optics.
[m,n] = size(inPSF);
if (n ~= m)
    error('Input must be a square matrix');
end
nLinearPixels = m;

[peakRow,peakCol] = FindMatPeak(inPSF);
radiusMat = MakeRadiusMat(nLinearPixels,nLinearPixels,peakCol,peakRow);
outPSF = zeros(nLinearPixels,nLinearPixels);
nBands = round(nLinearPixels/quantizationFactor);
radii = linspace(0,0.75*nLinearPixels,nBands);
for q = 1:length(radii)-1;
    index = find(radiusMat >= radii(q) & radiusMat < radii(q+1));
    if (~isempty(index))
        outPSF(index) = mean(inPSF(index));
    end
end
outPSF = outPSF/sum(outPSF(:));