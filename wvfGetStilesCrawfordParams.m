function sceParams = wvfGetStilesCrawfordParams(wls,source)
% sceParams = wvfGetStilesCrawfordParams(wls,source)
%
% Return a structure with the parameters to use
% for incorporating the Stiles-Crawford Effect
% into the Zernike optics calculations.
%
% Input
%   wls           -           Wavelengths (nm) over which to return rho
%   source        -           Data source
%
% Sources
%   'none'        -           Fill in rho with 0's leads to no SCE correction.
%   'berendshot'  -           Adjusted Berendshot et al. data for rho (default).
%
% Output
%   sceParams.wavelengths     SCE wavelengths.  This is just the passed wls.
%   sceParams.rho -           SCE rho as a function of wavelength.
%   sceParams.xo  -           SCE x center in mm
%   sceParams.yo              SCE y center in mm
%
% Berenschot data is taken from Berendshot et al. and then adjusted
% slighlty (subtracting .0045) to give rho=.041 at 550nm in agreement with
% Enoch and Lakshri's average foveal data in normals.  The xo, yo parameters
% are what Heidi had coded.
%
% If values for wavelengths outside of those over which data are specified
% are requested, the routine estimates by extending the last available value.
%
% Code provided by Heidi Hofer.
%
% Depends on Psychtoolbox splining routines.
%
% 8/21/11  dhb  Pulled into a separate routine.


wls = MakeItWls(wls);
switch (source)
    case 'berendshot'
        wls0 = (400:10:700)';
        rho0 = [0.0565 0.0585 0.0605 0.06 0.05875 0.05775 0.0565 0.0545 0.0525 0.051 0.04925 0.04675 0.044 0.041 0.04 0.041 0.041 0.0415 0.0425 0.04275 0.04325 0.045 0.047 0.048 0.0485 0.049 0.04975 0.05 0.04975 0.04925 0.049]';
        sceParams.xo=0.47;   % SCE center in mm
        sceParams.yo=0.00;   % SCE center in mm
    case 'none'
        wls0 = (400:10:700)';
        rho0 = zeros(size(wls0));
        sceParams.xo=0;   % SCE center in mm
        sceParams.yo=0;   % SCE center in mm
    otherwise
        error('Unsupported method passed');
        
end

% Spline to passed wavelenght sampling.
sceParams.wavelengths = wls;
sceParams.rho = SplineSrf(wls0,rho0,wls,1);
