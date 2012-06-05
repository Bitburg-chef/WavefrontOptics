function val = wvfGet(wvf,parm,varargin)
% val = wvfGet(wvf,parm,varargin)
%
% Get wavefront structure parameters and derived properties
%
% See also: wvfSet, wvfCreate, sceCreate, sceGet
%
% Wavefront properties are either stored as parameters or computed from those
% parameters. We generally store only unique values and  calculate all
% derived values.
%
%  A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used, where
%  unit specifies the spatial scale of the returned value:
%    length: 'm', 'cm', 'mm','um', 'nm'.
%    angle: 'deg', 'min', 'sec'
%
%  A leading '+ indicates that this is a get only parameter and may not be set.
%
% Parameters:
%
%  Bookkeeping
%   'name' - Name of this object
%   'type' - Type of this object, should always be 'wvf'
%
%  Zernike coefficients and measurement related
%   'zcoeffs' - Zernike coefficients
%   'measured pupil size' - Pupil size for wavefront aberration meaurements (mm,*)
%   'measured wl' - Wavefront aberration measurement wavelength (nm,*)
%   'measured optical axis' - Measured optical axis (deg)
%   'measured observer accommodation' - Observer accommodation at aberration measurement time (diopters)
%   'measured observer focus correction' - Focus correction added optically for observer at measurement time (diopters)
%
%  Spatial sampling parameters
%    'sample interval domain' - Which domain has sample interval held constant with wavelength ('psf', 'pupil')
%    'spatial samples' - Number of spatial samples (pixel) for pupil function and psf
%    'ref pupil plane size' - Size of sampled pupil plane at measurement wavelength (mm,*)
%    'ref pupil plane sample interval' - Pixel sample interval in pupil plane at measurement wavelength (mm,*)
%    'ref psf sample interval' - Sampling interval for psf at measurment wavelength (arcminute/pixel)
%  + 'pupil plane size' - Size of sampled pupil plane at any calculated wavelength(s) (mm)
%  + 'psf arcmin per sample' - Sampling interval for psf at any calculated wavelength(s) (min/pixel)
%  + 'psf angle per sample' - Sampling interval for psf at any calculated wavelength(s) (min,*/pixel)
%  + 'psf angular samples' - One-d slice of sampled angles for psf, centered on 0, for a single wavelength (min,*)
%
%  Spectral
%     'calc wavelengths' - Wavelengths to calculate over (nm,*)
%  +  'number calc wavelengths' - Number of wavelengths to calculate over
%
%     'weightspectrum'
%
% Pupil parameters
%     'calc pupil size'  - Pupil size for calculation (mm,*)
%     'calc optical axis' - Optical axis to compute for (deg)
%     'calc observer accommodation' - Observer accommodation at calculation time (diopters)
%     'calc observer focus correction' - Focus correction added optically for observer at calculation time (diopters)
%
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
DIDAGET = false;

%% Bookkeeping
switch parm
    case 'name'
        val = wvf.name;
        DIDAGET = true;
    case 'type'
        val = wvf.type;
        DIDAGET = true;
end

%% Zernike coefficients and related
switch parm
    case {'zcoeffs','zcoeff','zcoef'}
        % wvfGet(wvf,'zcoeffs',wllist)
        % wvfGet(wvf,'zcoeffs',4)
        if isempty(varargin),   val = wvf.zcoeffs;
        else                    val = wvf.zcoeffs(varargin{1});
        end
        DIDAGET = true;
        
    case {'measuredpupilsize', 'measuredpupil', 'measuredpupilmm', 'measuredpupildiameter'}
        % Pupil diameter in mm over for which wavefront expansion is valid
        % wvfGet(wvf,'measured pupil','mm')
        % wvfGet(wvf,'measured pupil')
        val = wvf.measpupilMM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'measuredwl', 'measuredwavelength'}
        % Measurement wavelength (nm)
        val = wvf.measWlNM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-9)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'measuredopticalaxis', 'measuredopticalaxisdeg'}
        % Measurement optical axis, degrees eccentric from fovea
        val = wvf.measOpticalAxisDeg;
        DIDAGET = true;
        
    case {'measuredobserveraccommodation', 'measuredobserveraccommodationdiopters'}
        % Observer accommodation, in diopters relative to relaxed state of eye
        val = wvf.measObserverAcommodationDiopters;
        DIDAGET = true;
        
    case {'measuredobserverfocuscorrection', 'measuredobserverfocuscorrectiondiopters'}
        % Focus correction added optically for observer at measurement time (diopters)
        val = wvf.measObserverAcommodationDiopters;
        DIDAGET = true;
end

%% Spatial sampling parameters
switch (parm)
    case {'sampleintervaldomain'}
        % What's held constant with calculated wavelength.
        % Choices are 'psf' and 'pupil'
        val = wvf.constantSampleIntervalDomain;
        DIDAGET = true;
        
    case {'spatialsamples', 'npixels', 'fieldsizepixels'}
        % Number of pixels that both pupil and psf planes are discretized
        % with.  This is a master value.
        val = wvf.nSpatialSamples;
        DIDAGET = true;
        
    case {'refpupilplanesize', 'refpupilplanesizemm', 'fieldsizemm'}
        % Total size of computed field in pupil plane.  This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.
        val = wvf.refSizeOfFieldMM;
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'refpupilplanesampleinterval', 'refpupilplanesampleintervalmm', 'fieldsamplesize', 'fieldsamplesizemmperpixel'}
        % Pixel sample interval of sample pupil field. This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.
        val = wvf.refSizeOfFieldMM/wvf.nSpatialSamples;
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'refpsfsampleinterval' 'refpsfarcminpersample', 'refpsfarcminperpixel'}
        % Arc minutes per pixel of the sampled psf at the measurement
        % wavelength.  This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.
        radiansPerPixel = wvfGet(wvf,'measured wl','mm')/wvfGet(wvf,'ref pupil plane size','mm');
        val = (180*60/3.1416)*radiansPerPixel;
        DIDAGET = true;
        
    case {'pupilplanesize', 'pupilplanesizemm'}
        % Total size of computed field in pupil plane, for calculated wavelengths(s).
        
        % Get wavelengths
        wavelengths = wvfGet(wvf,'calc wavelengths','nm');
        waveIdx = varargin{2};
        
        % Figure out what's being held constant with wavelength and act
        % appropriately.
        whichDomain = wvfGet(wvf,'sample interval domain');
        if (strcmp(whichDomain,'psf'))
            val = wvfGet(wvf,'ref pupil plane size','mm')*wavelengths(waveIdx)/wvfGet(wvf,'measured wl','nm');
        elseif (strcmp(whichDomain,'pupil'))
            val = wvfGet(wvf,'ref pupil plane size','mm')*ones(length(waveIdx),1);
        else
            error('Unknown sample interval domain ''%s''',whichDomain);
        end
        
        % Unit conversion
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'psfarcminpersample', 'psfarcminperpixel', 'arcminperpix'}
        % Arc minutes per pixel in psf domain, for calculated wavelength(s).
        
        % Get wavelengths
        wavelengths = wvfGet(wvf,'calc wavelengths','mm');
        waveIdx = varargin{1};
        
        % Figure out what's being held constant with wavelength and act
        % appropriately.
        whichDomain = wvfGet(wvf,'sample interval domain');
        if (strcmp(whichDomain,'psf'))
            val = wvfGet(wvf,'ref psf arcmin per pixel')*ones(length(waveIdx),1);
        elseif (strcmp(whichDomain,'pupil'))
            radiansPerPixel = wavelengths(waveIdx)/wvfGet(wvf,'ref pupil plane size','mm');
            val = (180*60/3.1416)*radiansPerPixel;
        else
            error('Unknown sample interval domain ''%s''',whichDomain);
        end
        DIDAGET = true;
        
    case {'psfanglepersample','angleperpixel','angperpix'}
        % Angular extent per pixel in the psf domain, for calculated wavelength(s).
        % wvfGet(wvf,'psf angle per sample',unit,waveIdx)
        %  unit = 'min' (default), 'deg', or 'sec'
        unit = varargin{1}; waveIdx = varargin{2};
        val = wvfGet(wvf,'psf arcmin per sample',waveIdx);
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
        DIDAGET = true;
        
    case {'psf angular samples','samplesangle','samplesarcmin','supportarcmin'}
        % Return one-d slice of sampled angles for psf, centered on 0, for a single wavelength
        % wvfGet(wvf,'psf angular samples',unit,waveIdx)
        %  unit = 'min' (default), 'deg', or 'sec'
        unit = varargin{1}; waveIdx = varargin{2};
        if (length(waveIdx) > 1)
            error('This only works for one wavelength at a time');
        end
        anglePerPix = wvfGet(wvf,'psf angle per sample',unit,waveIdx);
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'spatial samples');
        val = anglePerPix*((1:nPixels)-middleRow);
        DIDAGET = true;
end

%% Wavelength related
switch parm
    case {'calcwavelengths','wavelengths','wavelength','wls','wave'}
        % Wavelengths to compute on
        % wvfGet(wvf,'wave',unit,idx)
        % wvfGet(wvf,'wave','um',3)
        % May be a vector or single wavelength
        val = wvf.wls;
        
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*(1e-9)*ieUnitScaleFactor(unit);
        end
        
        % Select wavelength if indices were passed
        if length(varargin) > 1, val = val(varargin{2}); end
        DIDAGET = true;
        
    case {'numbercalcwavelengths','nwavelengths','nwave'}
        % Number of wavelengths to calculate at
        val = length(wvf.wls);
        DIDAGET = true;
        
%     case {'infocuswavelength','infocuswave','nominalfocuswl'}
%         % This should go away soon
%         val = wvfGet(wvf,'measured wl');
%         DIDAGET = true;
        
    case 'weightspectrum'
        val = wvf.weightingSpectrum;
        DIDAGET = true;    
end

%% Pupil parameters
switch parm
    case {'calcpupilsize', 'calculatedpupil'}
        % Pupil size to assume when computing pupil function and psf.  Must
        % be less than or equal to measured pupil size.
        %  wvfGet(wvf,'calculated pupil','mm')
        %  wvfGet(wvf,'calculated pupil','um')
        val = wvf.calcpupilMM;
        
        % Adjust units
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
    
      case {'calcopticalaxis'}
        % Specify optical axis at calculation time
        val = wvf.calcOpticalAxisDegrees;
        if (val ~= wvfGet(wvf,'measuredobserveraccommodation'))
            error('We do not currently know how to deal with values that differ from measurement time');
        end
        DIDAGET = true;
        
    case {'calcobserveraccommodation'}
        % Specify observer accommodation at calculation time
        val = wvf.calcObserverAccommodationDiopters;
        if (val ~= wvfGet(wvf,'measuredobserveraccommodation'))
            error('We do not currently know how to deal with values that differ from measurement time');
        end
        DIDAGET = true;
        
    case {'calcobserverfocuscorrection', 'defocusdiopters'}
        % Specify optical correction added to observer focus at calculation time
        val = wvf.calcObserverFocusCorrectionDiopters;
        DIDAGET = true;
        
    case {'pupilfunction','pupilfunc','pupfun'}
        % wvfGet(wvf,'pupilfunc',idx)  (idx <= nWave)
        %
        % The pupil function is derived from Zernicke coefficients in the
        % routine wvfComputePupilFunction If there are multiple
        % wavelengths, then this is a cell array of matrices. The sizes can
        % be a little different across wavelengths (see wvfComputePSF for
        % the relevant code). It has to do with scaling the pixel size to
        % be wavelength independent.  More explanation needed.
        
        % Can't do the get unless it has already been computed and is not stale.
        if (~isfield(wvf,'pupilfunc') || ~isfield(wvf,'PUPILFUNCTION_STALE') || wvf.PUPILFUNCTION_STALE)
            error('Must explicitly compute pupil function on wvf structure before getting it.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        % Return whole cell array of pupil functions over wavelength if
        % no argument passed.  If there is just one wavelength, we 
        % return the pupil function as a matrix, rather than as a cell
        % array with one entry.
        if isempty(varargin)
            if (length(wvf.pupilfunc) == 1)
                val = wvf.pupilfunc{1}; 
            else
                val = wvf.pupilfunc;
            end
        else
            idx = varargin{1}; nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.pupilfunc{idx};
            end
        end  
        DIDAGET = true;
        
        %     case {'defocusmicrons','defocusdistance'}
        %         % The defocus in distance rather than diopters
        %         % The default is microns.
        %         % wvfGet(wvfP,'defocus distance','mm');
        %         val = wvfGetDefocusFromWavelengthDifference(wvf);
        %         if ~isempty(varargin)
        %             % There is a different unit.  So, convert microns to meters and
        %             % then scale to new unit.
        %             val = (val/10^6)*ieUnitScaleFactor(varargin{1});
        %         end
        %         DIDAGET = true;
end

% Stiles Crawford Effect
switch parm
    case 'sceparams'
        if isfield(wvf,'sceParams'), val = wvf.sceParams; end
        DIDAGET = true;
        
    case 'scex0'
        if checkfields(wvf,'sceParams','xo'), val = wvf.sceParams.xo;
        else val = 0;
        end
        DIDAGET = true;
        
    case 'scey0'
        if checkfields(wvf,'sceParams','yo'), val = wvf.sceParams.yo;
        else val = 0;
        end
        DIDAGET = true;
        
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
        DIDAGET = true;
        
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
        DIDAGET = true;
        
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
        DIDAGET = true;
end
        
        
% Point and line spread data
switch parm
    case 'psf'
        % wvfGet(wvf,'psf',idx)  (idx <= nWave)
        
        % Force user to code to explicitly compute the psf if it isn't done.  Not ideal
        % but should be OK.
         if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute psf on wvf structure before getting it.  Use wvfComputePSF');
        end 
         
        % Return whole cell array of pupil functions over wavelength if
        % no argument passed.  If there is just one wavelength, we 
        % return the pupil function as a matrix, rather than as a cell
        % array with one entry.
        if isempty(varargin)
            if (length(wvf.psf) == 1)
                val = wvf.psf{1}; 
            else
                val = wvf.psf;
            end
        else
            idx = varargin{1}; nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.psf{idx};
            end
        end  
        DIDAGET = true;

% 'diffractionpsf' is broken right now, because of the way I redid the logic
% of the get on the psf.  Do we need it?  I think better to have
% a program that wants the diffraction limited PSF to build a
% a wfv object with zeros as the coefficients and simply operate
% on that.
%
% Breaking 'diffractionpsf' also breaks 'strehl'.  Hmmm.
%
%     case 'diffractionpsf'
%         % wvfGet(wvf,'diffraction psf',waveIdx);
%         % diffraction limited psf at wave(waveIdx)
%         %
%         waveIdx = varargin{1}; wave = wvfGet(wvf,'wave');
%         wvf = wvfSet(wvf,'wave',wave(waveIdx));
%         zcoeffs = zeros(65,1); zcoeffs(1) = 1;
%         wvf = wvfSet(wvf,'zcoeffs',zcoeffs);
%         wvf = wvfComputePSF(wvf);
%         val = wvfGet(wvf,'psf',1);
%         DIDAGET = true;
%
%     case 'strehl'
%         % wvfGet(wvf,'strehl',waveIdx);
%         % Strehl ratio.
%         % The strehl is the ratio of the peak of diff limited and the
%         % existing psf at that wavelength.
%         
%         % We could write this so that with no arguments we return all of
%         % the ratios across wavelengths.  For now, force a request for a
%         % wavelength index.
%         waveIdx = varargin{1};
%         psf = wvfGet(wvf,'psf',waveIdx);
%         dpsf = wvfGet(wvf,'diffraction psf',waveIdx);
%         val = max(psf(:))/max(dpsf(:));
%         
%         %         areaPixapod = wvfGet(wvf,'area pixapod',waveIdx);
%         %         val = max(psf(:))/areaPixapod^2;
%         %         % Old calculation was done in the compute pupil function routine.
%         % Now, we do it on the fly in here, for a wavelength
%         % strehl(wl) = max(max(psf{wl}))./(areapixapod(wl)^2);
%         DIDAGET = true;

    case 'psfcentered'
        % Centered so that peak is at middle position in coordinate grid
        val = psfCenter(wvfGet(wvf,'psf'));
        DIDAGET = true;
        
    case '1dpsf'
        % wvfGet(wvf,'1d psf',waveIdx,row)
        
        waveIdx = 1;
        whichRow = wvfGet(wvf,'middle row');
        if length(varargin) > 1, whichRow = varargin{2}; end
        if ~isempty(varargin),   waveIdx = varargin{1}; end
        
        psf = psfCenter(wvfGet(wvf,'psf',waveIdx));
        val = psf(whichRow,:);
        DIDAGET = true;
        
        
    case {'middlerow'}
        val = floor(wvfGet(wvf,'npixels')/2) + 1;
        DIDAGET = true;
        
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
        DIDAGET = true;
        
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
        DIDAGET = true;
        
    case {'distanceperpix','distperpix','distanceperpixel'}
        % Distance per pixel in specified unit ('mm')
        %   wvf(wvf,'distance per pixel','um');
        if isempty(varargin), unit = 'mm';
        else unit = varargin{1};
        end
        val = wvfGet(wvf,'field size',unit)/wvfGet(wvf,'npixels');
        DIDAGET = true;
        
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
        DIDAGET = true;
        
end

%% Catch the case where we don't know about the requested parameter
switch (parm)
    otherwise
        if (~DIDAGET)
            error('Unknown parameter %s\n',parm);
        end
end

return
