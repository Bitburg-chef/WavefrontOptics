function sceP = sceCreate(wave,source)
% sceP = sceCreate(wave,source)
%
% Return a structure with the Stiles-Crawford Effect parameters needed for the
% Zernike optics calculations.   
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
%   sceP.wavelengths     SCE wavelengths.  
%   sceP.rho -           SCE rho as a function of wavelength.
%   sceP.xo  -           SCE x center in mm
%   sceP.yo              SCE y center in mm
%
% Berendshot data is taken from Berendshot et al. and then adjusted
% slighlty (subtracting .0045) to give rho=.041 at 550nm in agreement with
% Enoch and Lakshri's average foveal data in normals.  The xo, yo parameters
% are what Heidi had coded.
%
% If values for wavelengths outside of those over which data are specified
% are requested, the routine estimates by extending the last available value.
%
% Code provided by Heidi Hofer.
%
% See also: sceGet, also depends on Psychtoolbox splining routines.
%
% Examples:
%    sceCreate
%    sceCreate(400:5:700,'berendshot')
%
% 8/21/11  dhb  Pulled into a separate routine.
%
% (c) WVF Toolbox Team 2011-2012

% TODO
%   There should also be an sceSet/Get
%

if ~exist('wave','var')   || isempty(wave), wave = (400:10:700)';  end
if ~exist('source','var') || isempty(source), source = 'none'; end
wave = wave(:);

% Let's include units.
switch (source)
    case 'berendshot'
        % Berendshot et al., wavelength parameters for rho0
        % These wavelengths are (400:10:700)';
        initWLS = [400,10,31];
        rho0 = [0.0565 0.0585 0.0605 0.06 0.05875 0.05775 0.0565 0.0545 0.0525 0.051 0.04925 0.04675 0.044 0.041 0.04 0.041 0.041 0.0415 0.0425 0.04275 0.04325 0.045 0.047 0.048 0.0485 0.049 0.04975 0.05 0.04975 0.04925 0.049]';
        sceP.xo=0.47;   % SCE center in mm
        sceP.yo=0.00;   % SCE center in mm
    case 'none'
        % Fill in with zeros
        initWLS = MakeItS(wave);
        rho0 = zeros(size(wave));
        sceP.xo=0;   % SCE center in mm
        sceP.yo=0;   % SCE center in mm
    otherwise
        error('Unsupported method %s\n',source);
        
end

% Spline initial wavelength sampling to request in wls
sceP.wavelengths = wave(:);

% Interpolate the rho values, if there are multiple wavelength samples.  This
% hasn't been tested much.
if length(rho0) > 1
    outWLS = MakeItS(wave);
    sceP.rho = SplineSrf(initWLS,rho0,outWLS,1);
else
    sceP.rho = rho0;
end

%% End
