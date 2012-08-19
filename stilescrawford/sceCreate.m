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
        % Berendshot et al. 2001, wavelength parameters for rho0.
        % These values were digitized from the bold curve in Figure
        % 2 of that paper, and then splined onto an evenly spaced
        % wavelength basis.
        initWls = [400 5 71];
        rho0 = [ ...
            0.0609    0.0609    0.0609    0.0608    0.0607    0.0606    0.0603    0.0601    0.0597    0.0593    0.0588    0.0582    0.0577    0.0571    0.0563    0.0554 ...
            0.0546    0.0537    0.0525    0.0515    0.0505    0.0494    0.0484    0.0473    0.0463    0.0455    0.0449    0.0444    0.0439    0.0437    0.0437    0.0437 ...    
            0.0441    0.0447    0.0454    0.0463    0.0472    0.0481    0.0487    0.0491    0.0496    0.0502    0.0508    0.0515    0.0522    0.0528    0.0534    0.0538 ...
            0.0541    0.0544    0.0547    0.0547    0.0547    0.0546    0.0543    0.0540    0.0536    0.0533    0.0529    0.0525    0.0522    0.0519    0.0515    0.0512 ...
            0.0510    0.0507    0.0506    0.0505    0.0505    0.0505    0.0506 ...
        ]';
    
        % These wavelengths are (400:10:700)' and were read off
        % one of the graphs in the Berendschot by eye by HH, and then
        % 0.045 was subtracted.  DB
        % is not sure which graph, since they don't obviously line
        % up with any of the data.  They may be the central tendency
        % by eye of the multiple curves in Figure 1, but since Figure 2
        % gives the mean of the data in Figure 1 (with some adjustment),
        % this isn't clearly the case.
        % initWLS = [400,10,31];
        % rho0 = [0.0565 0.0585 0.0605 0.06 0.05875 0.05775 0.0565 0.0545 0.0525 0.051 0.04925 0.04675 0.044 0.041 0.04 0.041 0.041 0.0415 0.0425 0.04275 0.04325 0.045 0.047 0.048 0.0485 0.049 0.04975 0.05 0.04975 0.04925 0.049]';
        
        
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
