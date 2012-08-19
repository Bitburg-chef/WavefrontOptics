function sceP = sceCreate(wave,source)
% sceP = sceCreate(wave,source)
%
% Return a structure with Stiles-Crawford Effect parameters.   
%
% Input
%   wls           -           Wavelengths (nm) over which to return rho
%   source        -           Data source
%
% Sources
%   'none'              -     Fill in rho with 0's leads to no SCE correction (default).
%   'berendschot data'  -     Adjusted Berendschot et al. (2001, JOSA A) mean psychophysical data for rho.
%   'berendschot model'  -    Adjusted Berendschot et al. (2001, JOSA A) model for rho.

% Output
%   sceP.wavelengths     SCE wavelengths.  
%   sceP.rho -           SCE peakedness rho as a function of wavelength (units: 1/mm^2)
%   sceP.xo  -           SCE x center in mm relative to pupil center
%   sceP.yo              SCE y center in mm relative to pupil center
%
% The 'berendschot' data/model is taken from Berendschot et al., 2001, "Wavelength
% dependence of the Stiles-Crawford ...", JOSA A, 18, 1445-1451 and then adjusted
% slighlty (subtracting .0045) to give rho = 0.041 at 550 nm in agreement with
% Enoch and Lakshminaranayan's average foveal data in normals (reference for this?).
% The model is the bold curve in Figure 2, which includes choroidal backscatter.
%
% The xo, yo parameters are from Applegate & Lakshminaranayan, "Parametric representation
% of Stiles-Crawford functions: normal variation of peak location and directionality",
% 1993, JOSA A, 1611-1623.
%
% If values for wavelengths outside of those over which data are specified
% are requested, the routine estimates by extending the last available value.
%
% Original code provided by Heidi Hofer.
%
% See also: sceGet, also depends on Psychtoolbox splining routines.
%
% Examples:
%    sceCreate
%    sceCreate(400:5:700,'berendschot')
%
% 8/21/11  dhb  Pulled into a separate routine.
%
% (c) WVF Toolbox Team 2011-2012

% TODO
%   There should also be an sceSet/Get
%   Could allow unit specification rather than default mm-based units.
%

if ~exist('wave','var')   || isempty(wave), wave = (400:10:700)';  end
if ~exist('source','var') || isempty(source), source = 'none'; end
wave = wave(:);

source = ieParamFormat(source);
switch (source)
    case 'berendschot_data'
        DIGITIZED = 0;
        if (DIGITIZED)
            % Berendschot et al. 2001, wavelength parameters for rho0.
            % These values were digitized from the mean psychophysical
            % data in Figure 2 of that paper, and then splined onto an evenly spaced
            % wavelength basis.  The data are not very smooth. Heidi's
            % values taken by eye from the graph smooth it out some and
            % are probably more realistic.
            %
            % Units of peakedness in the data file are 1/mm^2.
            dataFilename = fullfile(wvfRootPath,'data','berendschotEtAl2001_Figure2Data.txt');
            rawData = ReadStructsFromText(dataFilename);
            initWls = [400 5 71];
            rho0 = interp1([rawData.Wavelength]',[rawData.Peakedness]',SToWls(initWls),'linear');
            index = find(SToWls(initWls) == 550);
            if (isempty(index))
                error('Oops.  Need a splined value at exactly 550 nm');
            end
            rho0 = rho0 - rho0(index) + 0.041;
            
        else 
            % These were the psychophysical data read off Figure 2 by eye by HH, with 0.045
            % subtracted to produce 0.041 at 550 nm.  I think they are probably preferable to the
            % digitized version above because they smooth the data.  We could figure out how to smooth
            % the psychophysical data from the digitized values, but not today.
            initWls = [400,10,31];
            rho0 = [0.0565 0.0585 0.0605 0.06 0.05875 0.05775 0.0565 0.0545 0.0525 0.051 0.04925 0.04675 0.044 0.041 0.04 0.041 0.041 0.0415 0.0425 0.04275 0.04325 0.045 0.047 0.048 0.0485 0.049 0.04975 0.05 0.04975 0.04925 0.049]';
        end

    case 'berendschot_model'
        % Berendschot et al. 2001, wavelength parameters for rho0.
        % These values were digitized from the bold curve in Figure
        % 2 of that paper, and then splined onto an evenly spaced
        % wavelength basis.
        %
        % Units of peakedness in the data file are 1/mm^2.
        dataFilename = fullfile(wvfRootPath,'data','berendschotEtAl2001_Figure2BoldSmooth.txt');
        rawData = ReadStructsFromText(dataFilename);
        initWls = [400 5 71];
        rho0 = interp1([rawData.Wavelength]',[rawData.Peakedness]',SToWls(initWls),'linear');
        index = find(SToWls(initWls) == 550);
        if (isempty(index))
            error('Oops.  Need a splined value at exactly 550 nm');
        end
        rho0 = rho0 - rho0(index) + 0.041;
              
    case 'none'
        % Fill in with zeros
        initWls = MakeItS(wave);
        rho0 = zeros(size(wave));

    otherwise
        error('Unsupported method %s\n',source);
        
end

sceP.xo=0.47;   % SCE center in mm
sceP.yo=0.00;   % SCE center in mm
        
        
% Spline initial wavelength sampling to request in wls,
% and extend with values at min/max wavelength if requested
sceP.wavelengths = wave(:);
sceP.rho = interp1(SToWls(initWls),rho0,sceP.wavelengths,'linear');
index = find(sceP.wavelengths < initWls(1));
if (~isempty(index))
    sceP.rho(index) = rho0(1);
end
index = find(sceP.wavelengths < initWls(end));
if (~isempty(index))
    sceP.rho(index) = rho0(end);
end


%% End
