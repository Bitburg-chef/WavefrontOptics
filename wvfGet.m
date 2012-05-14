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
%     'field sample size mm'
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
% History:
%   4/29/12  dhb  Allow wls as a synonym for wavelength, because that was
%                 my first guess given the name of the field.
%
% (c) Wavefront Toolbox Team 2011, 2012

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
    case {'wave','wavelength','wavelengths','wls'}
        % wvfGet(wvf,'wave',unit)
        % wvfGet(wvf,'wave','um')
        % May be a vector or single wavelength
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
    case {'infocuswavelength','infocuswave','nominalfocuswl'}
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
        % wvfGet(wvf,'measured pupil','mm')
        % wvfGet(wvf,'measured pupil')
        val = wvf.measpupilMM;               % Default pupil diameter?
        % Scale for unit
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        
        % Focus parameters
    case {'zcoeffs','zcoeff','zcoef'}
        % wvfGet(wvf,'zcoef',list)
        % wvfGet(wvf,'zcoef',4)
        if isempty(varargin),   val = wvf.zcoeffs;
        else                    val = wvf.zcoeffs(varargin{1});
        end
    case {'pupilfunction','pupilfunc','pupfun'}
        % wvfGet(wvf,'pupilfunc',idx)  (idx <= nWave)
        %
        % The pupil function is derived from Zernicke coefficients in the
        % routine wvfComputePupilFunction If there are multiple
        % wavelengths, then this is a cell array of matrices. The sizes can
        % be a little different across wavelengths (see wvfComputePSF for
        % the relevant code). It has to do with scaling the pixel size to
        % be wavelength independent.  More explanation needed.
        if isempty(varargin)
            % This is the whole cell array (if there are multiple) or just
            % the single matrix if there is only one wavelength.
            if isfield(wvf,'pupilfunc'), val = wvf.pupilfunc; end
        else
            idx = varargin{1}; nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.pupilfunc{idx};
            end
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
        else warning('WVFGET:scefract','No sceFrac field');
        end
        
        % Point and line spread data
    case 'psf'
        % wvfGet(wvf,'psf',idx)  (idx <= nWave)
        % The point spread function is calculated from the pupilfunction. I
        % almost think we should not store it, but always calculate it.
        % Computers are fast enough to deal with the fft, and the whole
        % computation is just a few lines long.
        if isempty(varargin)
            % This is the whole cell array (if there are multiple) or just
            % the single matrix if there is only one wavelength.
            if isfield(wvf,'psf'), val = wvf.psf; end
        else
            idx = varargin{1}; nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.psf{idx};
            end
        end
    case 'psfcentered'
        % Centered so that peak is at middle position in coordinate grid
        val = psfCenter(wvfGet(wvf,'psf'));
    case '1dpsf'
        % wvfGet(wvf,'1d psf',waveIdx,row)
        
        waveIdx = 1;
        whichRow = floor(wvfGet(wvf,'npixels')/2) + 1;
        if length(varargin) > 1, whichRow = varargin{2}; end
        if ~isempty(varargin),   waveIdx = varargin{1}; end
        
        psf = psfCenter(wvfGet(wvf,'psf',waveIdx));
        val = psf(whichRow,:);
    case 'strehl'
        % Strehl ratio.  Not sure when it is calculated
        % It seems to be wrong mostly.  We should fix it.
        % The strehl is the ratio of the peak of diff limited and the
        % existing.
        if isfield(wvf,'strehl'),  val = wvf.strehl;
        else                       disp('No strehl parameter present');
        end
        
        % Spatial and angular support
    case {'fieldsizepixels','npixels'}
        val = wvf.sizeOfFieldMM/wvf.fieldSampleSizeMMperPixel;
        if   val ~= round(val)
            warning('WVFGET:npixels','npixels not an integer.');
        else val = round(val);
        end
        
    case {'fieldsamplesize','fieldsamplesizemm'}
        val = wvf.fieldSampleSizeMMperPixel;
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
        
        % These pixel related measures used to be computed in wvfComputePSF
        % and wvfComputePupilFunction.
    case {'areapix'}
        % Not sure about the physical significance of this
        % If we know the area of each pixel, I suppose we can calculate the
        % area covered by the pupil function.
        if isempty(varargin)
            nWave = wvfGet(wvf,'n wave');
            val = zeros(nWave,1);
            for ii = 1:nWave
                val(ii) = numel(wvfGet(wvf,'pupil function',ii));
            end
        else
            val = numel(wvfGet(wvf,'pupil function',varargin{1}));
        end
        
    case {'areapixapod'}
        % Not sure about the physical significance of this
        % Something like the area underneath the absolute value of the
        % pupil function.
        if isempty(varargin)
            nWave = wvfGet(wvf,'n wave');
            val = zeros(nWave,1);
            for ii = 1:nWave
                val(ii) = sum(sum(abs(wvfGet(wvf,'pupil function',ii))));
            end
        else 
            val = sum(sum(abs(wvfGet(wvf,'pupil function',varargin{1}))));
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
        % Get the samples in degrees
        val = wvfGet(wvf,'samples angle','deg');
        
        % Convert to meters and then to selected spatial scale
        val = val*(300*10^-6);  % Sample in meters assuming 300 um / deg
        val = val*ieUnitScaleFactor(unit);

        % Old code.  This used the distance in the pupil plane, which is
        % wrong (I think).  I think that the angle calculation is probably
        % correct.  We should use the fact that we know that 1 deg in the
        % human eye is 300 um, and then get the better calculation (based
        % on numerical aperture?) from HH or DHB or someone.
        %         distPerPix = wvfGet(wvf,'distperpix',unit);
        %         middleRow  = wvfGet(wvf,'middle row');
        %         nPixels    = wvfGet(wvf,'npixels');
        %         val = distPerPix*((1:nPixels)-middleRow);
        
    otherwise
        error('Unknown parameter %s\n',parm);
end

return
