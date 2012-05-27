function  val = wvfGet(wvf,parm,varargin)
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
%Parameter list
%
% General
%     'name' - Name of this structure
%     'type' - Always 'wvf'
%
% Spectral
%     'wavelength' - wavelength samples     (wvfGet(wvfP,'wave',unit,idx))
%     'nwave'      - number of wavelengths  (wvfGet(wvfP,'n wave'))
%     'infocus wavelength'
%     'weightspectrum'
%
% Pupil parameters
%     'calculated pupil'  - Pupil size for calculation (mm)
%     'measured pupil'    - Pupil size when measured   (mm)
%     'field sample size' - Aperture .... mm
%     'field size pixels'
%     'field size'        - Aperture size (mm)
%
% Focus parameters
%     'zcoef'              - Zernicke polynomial coefficients (n=65)
%     'defocusdiopters'
%     'defocus distance'   -         *microns
%     'weight spectrum'

% Stiles Crawford Effect
%     'sce params' - The whole structure
%     'sce x0'
%     'sce y0'
%     'sce rho'
%     'sce wavelengths'*
%
% Pointspread function
%     'psf'            - Point spread function
%     'diffractionpsf' - Diffraction limited psf for these parameters
%     'psf centered'   - 
%     '1d psf'
%     'spatialsupport'
%     'middle row'
%     'samples space'  - Includes key parameter umPerDeg, needs thought
%     'samples angle'
%
%     'arcminperpix'
%     'angleperpixel'
%     'strehl'     - Ratio of peak of diffraction limited to actual
%
%     'areapixapod' - In need of repair to speed up strehl
%     'areapix' - In need of repair to speed up strehl
%
% Examples:
% * Compute diffraction limited psf
%   wvfP = wvfCreate; 
%   wvfP = wvfComputePSF(wvfP); 
%   vcNewGraphWin; wvfPlot(wvfP,'image psf','um',1,20)
%   psf = wvfGet(wvfP,'diffraction psf',1); vcNewGraphWin; mesh(psf)
%
% * Strehl is ratio of diffraction and current
%   wvfP = wvfComputePSF(wvfP); wvfGet(wvfP,'strehl',1)
%
% * Blur and recompute.  4th coefficient is defocus
%   z = wvfGet(wvfP,'zcoeffs');z(4) = 0.3; wvfP = wvfSet(wvfP,'zcoeffs',z); 
%   wvfP = wvfComputePSF(wvfP); wvfGet(wvfP,'strehl',1)
%
% See also: wvfComputePupilFunction, wvfLoadHuman
%
% (c) Wavefront Toolbox Team 2011, 2012
%
% History
%   5/22/12 dhb      Improve comment about how arcminperpix is obtained.

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
        % wvfGet(wvf,'wave',unit,idx)
        % wvfGet(wvf,'wave','um',3)
        % May be a vector or single wavelength
        val = wvf.wls;
        
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*(1e-9)*ieUnitScaleFactor(unit);
        end
        % A selected wavelength
        if length(varargin) > 1, val = val(varargin{2}); end
        
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
        if ~isempty(varargin)
            % There is a different unit.  So, convert microns to meters and
            % then scale to new unit.
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
        
    case {'scefrac','scefraction','stilescrawfordeffectfraction'}
        % wvfGet(wvf,'sce fraction',waveIdx)
        % The variable sceFrac tells you how much light
        % is lost if you correct for the SCE.
        waveIdx = varargin{1};
        val = wvfGet(wvf,'area pixapod',waveIdx) / wvfGet(wvf,'areapix',waveIdx);
        
        % Note what this is or when it is calculated
        %         if checkfields(wvf,'sceFrac'), val = wvf.sceFrac;
        %         else warning('WVFGET:scefract','No sceFrac field');
        %         end
        
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
            else
                if checkfields(wvf,'psf'), val = wvf.psf{idx};
                else disp('No psf.  Use wvfComputePSF');
                end
            end
        end
    case 'diffractionpsf'
        % wvfGet(wvf,'diffraction psf',waveIdx);
        % diffraction limited psf at wave(waveIdx)
        %
        waveIdx = varargin{1}; wave = wvfGet(wvf,'wave'); 
        wvf = wvfSet(wvf,'wave',wave(waveIdx));
        zcoeffs = zeros(65,1); zcoeffs(1) = 1;
        wvf = wvfSet(wvf,'zcoeffs',zcoeffs); 
        wvf = wvfComputePSF(wvf);
        val = wvfGet(wvf,'psf',1);
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
        % wvfGet(wvf,'strehl',waveIdx);
        % Strehl ratio.  
        % The strehl is the ratio of the peak of diff limited and the
        % existing psf at that wavelength.
        
        % We could write this so that with no arguments we return all of
        % the ratios across wavelengths.  For now, force a request for a
        % wavelength index.
        waveIdx = varargin{1};
        psf = wvfGet(wvf,'psf',waveIdx);
        dpsf = wvfGet(wvf,'diffraction psf',waveIdx);
        val = max(psf(:))/max(dpsf(:));
        
        %         areaPixapod = wvfGet(wvf,'area pixapod',waveIdx);
        %         val = max(psf(:))/areaPixapod^2;
        %         % Old calculation was done in the compute pupil function routine.
        % Now, we do it on the fly in here, for a wavelength
        % strehl(wl) = max(max(psf{wl}))./(areapixapod(wl)^2);
        
        % Spatial and angular support
    case {'fieldsizepixels','npixels'}
        val = wvf.sizeOfFieldMM/wvf.fieldSampleSizeMMperPixel;
        if   val ~= round(val)
            warning('WVFGET:npixels','npixels not an integer.');
        else val = round(val);
        end
        
    case {'fieldsamplesize','fieldsamplesizemm'}
        % wvfGet('field sample size','mm',waveIdx)
        %
        % This quantity is scaled to make the field size the same,
        % independent of wavelength. It is always normalized for the
        % setScaleWl of 550nm.  This should be explained in terms of
        % physics better.
        % The value stored in field samplesizeMMperPixel is for 550nm, I
        % guess. 
        % OLD CODE:  val = wvf.fieldSampleSizeMMperPixel;
        
        % For now, force the long call.  Otherwise an error.
        unit = varargin{1};  waveIdx = varargin{2};
        setScaleWl = 550; wave = wvfGet(wvf,'wave');
        %  if ~isempty(varargin), unit = varargin{1}; end
        %  if length(varargin) > 1, waveIdx = varargin{2}; end
        val = wvf.fieldSampleSizeMMperPixel*(wave(waveIdx)/setScaleWl);
        
        if ~isempty(unit)
            % Convert from mm to meters and then scale to unit
            val = (val/1000)*ieUnitScaleFactor(unit);
        end
    case {'angleperpixel','angperpix'}
        % wvfGet(wvf,'angle per pixel',unit,waveIdx)
        % Angle per pixel in various angle units
        %  unit = 'min', 'deg', or 'sec'
        unit = varargin{1}; waveIdx = varargin{2};
        % The angle units are always in minutes
        val = wvfGet(wvf,'arcminperpix',waveIdx);

        if ~isempty(unit)
            unit = lower(unit);
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
        % wvfGet(wvf,'arcmin per pix',waveIdx)
        % I think the pixels here are in the pupil plane.
        wave = wvfGet(wvf,'wave');
        waveIdx = varargin{1};
        thisWave = wave(waveIdx);   % Wavelength in nanometers
        pupilPlaneFieldMM = wvfGet(wvf,'field size','mm',waveIdx);
        
        % When we convert between the pupil function and the PSF,
        % we use the fft.  Thus the size of the image in pixels 
        % is the same for the sampled pupil function and the sampled
        % psf.
        %
        % The number of arc minutes per pixel in the sampled PSF is
        % related to the number of mm per pixel for hte pupil function,
        % with the relation depending on the wavelength.  The fundamental
        % formula in the pupil plane is that the pixel sampling interval
        % in cycles/radian is:
        %
        %   pupilPlaneCyclesRadianPerPix = pupilPlaneField/[lambda*npixels]
        %
        % where npixels is the number of linear pixels and lambda is the
        % wavelength. This formula may be found as Eq 10 of Ravikumar et al.
        % (2008), "Calculation of retinal image quality for polychromatic light,"
        % JOSA A, 25, 2395-2407, at least if we think their quantity d is the
        % size of the pupil plane field being sampled.
        %
        % If we now remember how units convert when we do the fft, we obtain
        % that the number of radians in the PSF image is the inverse of the
        % sampling interval:
        %
        %   radiansInPsfImage = [lambda*npixels]/pupilPlaneField
        %
        % which then gives us the number of radiansPerPixel in the
        % PSF image as
        %
        %   radiansPerPixel = lambda/pupilPlaneField
        % 
        % The formula below implements this, with a conversion
        % from radians to minutes with factor (180*60/3.1416)
        % and converts wavelength to mm from nm with factor (.001*.001)
        %
        % DHB, 5/22/12, based on earler comments that were here.  Someone
        % else might take a look at the paper referenced above and the logic
        % of this comment and check that it all seems right.  Did I think
        % through the fft unit conversion correctly?  And, there must be
        % a more fundamental reference than the paper above, and for which
        % one wouldn't have to guess quite as much about what is meant.
        val = (180*60/3.1416)*thisWave*(.001*.001)/pupilPlaneFieldMM;
        
    case {'samplesangle','samplesarcmin','supportarcmin'}
        %    wvfGet(wvf,'samples angle','min',waveIdx)
        %
        % Sample support in angle ('min' default), centered on 0
        %  angle can also be 'deg' or 'sec'
        % This is apparently wavelength dependent.  We should use the same
        % method here as we use in wvfComputePSF.  
        unit = varargin{1}; waveIdx = varargin{2};
        anglePerPix = wvfGet(wvf,'angleperpixel',unit,waveIdx);
        
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'npixels');
        val = anglePerPix*((1:nPixels)-middleRow);
        
    case {'middlerow'}
        val = floor(wvfGet(wvf,'npixels')/2) + 1;
        
    % Used to be assigned, but now computed.
    case {'fieldsize','fieldsizemm','fieldsizespace'}
        % wvfGet(wvf,'field size','mm',waveIdx)
        % The field size in the pupil plane in mm
        % The size of a pixel * the number of pixels is the field size.
        % We used to store all three, but that could lead to inconsistency.
        unit = varargin{1}; waveIdx = varargin{2};
        npixels = wvfGet(wvf,'npixels');
        val = wvfGet(wvf,'field sample size',unit,waveIdx)*npixels;
        % OLD: val = wvf.sizeOfFieldMM;  % This parameter should go away.
        
        % The two measures below used to be computed and assigned in
        % wvfComputePSF and wvfComputePupilFunction.  They were win wvfSet.
        % But now we just compute them on the fly, here. - BW
    case {'areapix'}
        % Not sure about the physical significance of this If we know the
        % area of each pixel, we can use this to calculate the area covered
        % by the pupil function. It is computed for the first time in
        % wvfComputePupilFunction as numel(pupilfunc))
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
        % Not sure about the physical significance of this Something like
        % the area underneath the absolute value of the pupil function It
        % is a vector the same length as wavelength. It is computed for the
        % first time in wvfComputePupilFunction sum(sum(abs(pupilfunc)))
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
        % wvfGet(wvf,'samples space','um')
        % Spatial support in samples, centered on 0
        % Unit and wavelength must be specified
        
        % This parameter matters for the OTF and PSF quite a bit.
        umPerDeg = (330*10^-6);
        
        unit = varargin{1}; waveIdx = varargin{2};
        % Get the samples in degrees
        val = wvfGet(wvf,'samples angle','deg',waveIdx);
        
        % Convert to meters and then to selected spatial scale
        val = val*umPerDeg;  % Sample in meters assuming 300 um / deg
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
