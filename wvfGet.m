function val = wvfGet(wvf,parm,varargin)
%Get wavefront structure parameters and derived properties
%
%     val = wvfGet(wvf,parm,varargin)
%
% Wavefront properties are either stored as parameters or computed from those
% parameters. We generally store only unique values and  calculate all
%  derived values.
%
%  A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used, where
%  unit specifies the spatial scale of the returned value:  'm', 'cm', 'mm',
%  'um', 'nm'.  Default is always meters ('m').  (REALLY?  NOT YET).
%
% Examples:
%
%
%Parameters
% General
%     'name' - Name of this wavefront parameter structure
%     'type' - Always 'wvf'
%
% Spectral 
%     'wavelength' - wavelength samples
%     'nwave' - Number of wavelength samples
%     'infocus wavelength'
%
% Pupil parameters
%     'calculated pupil'
%     'measured pupil'
%
% Field?
%     'field size pixels'
%     'field size mm'
%
% Focus parameters
%     'zcoef'
%     'defocusdiopters'
%     'defocus distance' *microns
%     'weight spectrum'

% Stiles Crawford Effect
%     'sce params' - The whole structure
%     'sce x0'
%     'sce y0'
%     'sce rho'
%     'sce wavelengths'*
%
% Pointspread function
%     'psf'
%     'psf centered'
%     '1d psf'
%     'strehl'
%
% (c) Wavefront Toolbox Team 2011

if ~exist('parm','var') || isempty(parm), error('Parameter must be defined.'); end

% Default is empty when the parameter is not yet defined.
val = [];

parm = ieParamFormat(parm);

switch parm
    case 'name'
        val = wvf.name;
    case 'type'
        val = wvf.type;
   
        % Spectral matters
    case {'wave','wavelength','wavelengths'}
        % wvfGet(wvf,'wave',unit)
        % wvfGet(wvf,'wave','um')
        val = wvf.wls;
        
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*(1e-9)*ieUnitScaleFactor(unit);
        end

    % Wavelength related
    case 'weightspectrum'
        val = wvf.weightingSpectrum;         % Defocus
    case 'nwave'
        val = length(wvf.wls);
    case {'infocuswavelength','infocuswave'}
        val = wvf.nominalFocusWl;            % In focus wavelength (nm)
        
        % Pupil parameters
    case 'calculatedpupil'
        % Measured describes original data.  Calculated describes what we
        % are using in the simulation.  Default in mm.
        %  wvfGet(wvf,'calculated pupil','mm')
        %  wvfGet(wvf,'calculated pupil','um')
        val = wvf.calcpupilMM;               % Default pupil? diameter?
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
    case 'measuredpupil'
        % Measurements - maximum value in mm
        % Default is in millimeters
        % wvfGet(wvf,'measured pupil','m')
        % wvfGet(wvf,'measured pupil')
        val = wvf.measpupilMM;               % Default pupil diameter?
        % Scale for unit
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
             
        % Focus parameters
    case {'zcoeffs'}
        % wvfGet(wvf,'zcoef',list)
        % wvfGet(wvf,'zcoef',4)
        if isempty(varargin)
            val = wvf.zcoeffs;
        else
            val = wvf.zcoeffs(varargin{1});
        end
    case 'defocusdiopters'
        val = wvf.defocusDiopters;           % Defocus
    case {'defocusmicrons','defocusdistance'}
        % The defocus in distance rather than diopters
        % The default is microns.
        % wvfGet(wvfP,'defocus distance','mm');
        val = wvfGetDefocusFromWavelengthDifference(wvf);
        if isempty(varargin), return
        else
            % Convert to meters and then scale
            val = (val/10^6)*ieUnitScaleFactor(varargin{1});
        end
              
        % Stiles Crawford Effect
    case 'sceparams'
        if isfield(wvf,'sceParams'), val = wvf.sceParams; end
    case 'scex0'
        if checkfields(wvf,'sceParams','xo'), val = wvf.sceParams.xo;
        else val = 0; 
        end
    case 'scey0'
        if checkfields(wvf,'sceParams','yo'), val = wvf.sceParams.yo;
        else val = 0;
        end
    case {'scewavelength','scewavelengths','scewave'}
        % This returns the wvf wavelength list if there isn't a sceParams
        % structure.  Might be OK.
        % wvfGet(wvf,'sce wavelengths',unit)
        if checkfields(wvf,'sceParams','wavelengths'), val = wvf.sceParams.wavelengths;
        else val = wvf.wls;
        end
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*10e-9*ieUnitScaleFactor(unit);
        end
    case 'scerho'
        % Get rho value for a particular wavelength
        %  wvfGet(wvf,'rho',waveList)
        if checkfields(wvf,'sceParams','rho'), val = wvf.sceParams.rho;
        else val = zeros(wvfGet(wvf,'nWave'),1);
        end
        
        % Return rho values for selected wavelengths
        if ~isempty(varargin)
            wave = wvfGet(wvf,'sce wave');  % The waves for rho
            waveList = varargin{1};
            index = find(wave == waveList);
            if ~isempty(index), val = val(index);
            else error('Passed wavelength not contained in sceParams');
            end
        end
    
    case 'scefrac'
        % Note what this is or when it is calculated
        if checkfields(wvf,'sceFrac'), val = wvf.sceFrac; 
        else warning('No sceFrac field');
        end
        
       % Point and line spread data
    case 'psf'
        val = wvf.psf;
    case 'psfcentered'
        % Centered so that peak is at middle position in coordinate grid 
        val = psfCenter(wvfGet(wvf,'psf'));
    case '1dpsf'
        % wvfGet(wvf,'1d psf',row)
        psf = psfCenter(wvfGet(wvf,'psf'));

        if isempty(varargin)
            whichRow = floor(wvfGet(wvf,'npixels')/2) + 1;
        else
            whichRow = varargin{1};
        end
        val = psf(whichRow,:);
     case 'strehl'
        % Strehl ratio.  Not sure when it is calculated
        if isfield(wvf,'strehl'),  val = wvf.strehl;
        else                       disp('No strehl parameter present');
        end
   
        % Spatial and angular support
    case {'fieldsizepixels','npixels'}
        % In pixels?  No units?  Why not a distance or an angle or
        % something?
        val = wvf.sizeOfFieldPixels;
    case {'angleperpixel','angperpix'}
        % wvfGet(wvf,'angle per pixel',unit)
        %  unit = 'min', 'deg', or 'sec'
        val = wvfGet(wvf,'arcminperpix');
        if ~isempty(varargin)
            unit = lower(varargin{1});
            switch unit
                case 'deg'
                    val = val/60;
                case 'sec'
                    val = val*60;
                case 'min'
                    % Default
                otherwise
                    error('unknown angle unit %s\n',unit);
            end
        end
        
    case 'arcminperpix'
        % Visual angle (arc minutes) per pixel.
        if isempty(varargin),  wave = 550;
        else                   wave = varargin{1};
        end
        fieldMM = wvfGet(wvf,'field size mm');
        % Wavelength in nanometers
        val = (180*60/3.1416)*wave*.001*.001/fieldMM;
    case {'samplesangle','samplesarcmin','supportarcmin'}
        % Sample support in angle ('min' default), centered on 0
        %    wvfGet(wvf,'samples angle','min')
        %  angle can also be 'deg' or 'sec'
        
        unit = 'min';
        if ~isempty(varargin),unit = varargin{1}; end
        anglePerPix = wvfGet(wvf,'angleperpixel',unit);
        
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'npixels');
        val = anglePerPix*((1:nPixels)-middleRow);
    case {'middlerow'}
        val = floor(wvfGet(wvf,'npixels')/2) + 1;    
    case {'fieldsizemm','fieldsizespace','fieldsize'}
        val = wvf.sizeOfFieldMM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val/1000)*ieUnitScaleFactor(varargin{1});
        end
    case {'distanceperpix','distperpix','distanceperpixel'}
        % Distance per pixel in specified unit ('mm')
        %   wvf(wvf,'distance per pixel','um');
        if isempty(varargin), unit = 'mm';
        else unit = varargin{1};
        end
        val = wvfGet(wvf,'field size',unit)/wvfGet(wvf,'npixels');
    case {'samplesspace','supportspace','spatialsupport'}
        % Spatial support in samples, centered on 0
        % Unit can be specified
        %    wvfGet(wvf,'samples space','um')
        if isempty(varargin), unit = 'mm';
        else unit = varargin{1};
        end
        distPerPix = wvfGet(wvf,'distperpix',unit);
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'npixels');
        val = distPerPix*((1:nPixels)-middleRow);
        
    otherwise
        error('Unknown parameter %s\n',parm);
end

return
