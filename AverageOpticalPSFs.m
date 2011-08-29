function [averagePSF] = AverageOpticalPSFs(inputPSFs)
% [averagePSF] = AverageOpticalPSFs(inputPSFs)
%
% Produces an average of passed optical point spread functions.
% Input is a cell array of PSFs.
% Average in sf domain and transform back.
%
% 7/19/07  dhb  Wrote it.

nInputs = length(inputPSFs);
for i = 1:nInputs
    inputOTFs{i} = psf2otf(inputPSFs{i});
end
averageOTF = zeros(size(inputOTFs{1}));
for i = 1:nInputs
    averageOTF = averageOTF + inputOTFs{i};
end
averageOTF = averageOTF/nInputs;
averagePSF = otf2psf(averageOTF);

