    
% Set values in millimeters
wvfParams0 = wvfCreate('measured pupil',6,'calculated pupil',3);

% Sample data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = load(sDataFile);

whichSubjects = 1:2;
theZernikeCoeffs = theZernikeCoeffs(:,whichSubjects);
nSubjects = size(theZernikeCoeffs,2);
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Stiles Crawford
DOSCE = 0;
theWavelength = 550;
if (DOSCE), wvfParams0.sceParams = sceCreate(theWavelength,'berendshot');
else        wvfParams0.sceParams = sceCreate(theWavelength,'none');
end

wvfParams0 = wvfComputePupilFunction(wvfParams0);
pupilfunc = wvfGet(wvfParams0,'pupil function');
vcNewGraphWin; imagesc(pupilfunc)

for ii=whichSubjects
    wvfParams = wvfSet(wvfParams0,'zcoeffs',theZernikeCoeffs(:,ii));
    wvfParams = wvfComputePupilFunction(wvfParams);
    vcNewGraphWin; mesh(angle(pupilfunc)); colormap(gray)
end



