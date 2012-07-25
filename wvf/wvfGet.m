function val = wvfGet(wvf,parm,varargin)
% val = wvfGet(wvf,parm,varargin)
%
% Get wavefront structure parameters and derived properties
%
% See also: wvfSet, wvfCreate, wvfComputePupilFunction, wvfComputePSF, sceCreate, sceGet
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
%    'number spatial samples' - Number of spatial samples (pixel) for pupil function and psf
%    'ref pupil plane size' - Size of sampled pupil plane at measurement wavelength (mm,*)
%    'ref pupil plane sample interval' - Pixel sample interval in pupil plane at measurement wavelength (mm,*)
%    'ref psf sample interval' - Sampling interval for psf at measurment wavelength (arcminute/pixel)
%  + 'pupil plane size' - Size of sampled pupil plane at any calculated wavelength(s) (mm)
%  + 'psf arcmin per sample' - Sampling interval for psf at any calculated wavelength(s) (min/pixel)
%  + 'psf angle per sample' - Sampling interval for psf at any calculated wavelength(s) (min,*/pixel)
%  + 'psf angular samples' - One-d slice of sampled angles for psf, centered on 0, for a single wavelength (min,*)
%  + 'psf spatial samples' - One-d slice of sampled psf in spatial units, centered on 0 for a single wavelength (*)
%  + 'pupil spatial samples' - One-d slice of sampled pupil function in spatial units, centered on 0 for a single wavelength (*)
%  + 'middle row' - The middle row of sampled functions
%
%  Calculation parameters
%     'calc pupil size'  - Pupil size for calculation (mm,*)
%     'calc optical axis' - Optical axis to compute for (deg)
%     'calc observer accommodation' - Observer accommodation at calculation time (diopters)
%     'calc observer focus correction' - Focus correction added optically for observer at calculation time (diopters)
%     'calc wavelengths' - Wavelengths to calculate over (nm,*)
%     'calc cone psf info' - Structure with cone sensitivities and weighting spectrum for computing cone psfs.
%  +  'number calc wavelengths' - Number of wavelengths to calculate over
%
% Stiles Crawford Effect
%     'sce params' - The whole structure
%     'sce x0'
%     'sce y0'
%     'sce rho'
%     'sce wavelengths'*
%  +  'sce fraction' - How much light is effectively lost by cones because of sce
%  +  'areapix' - Used in computation of sce fraction
%  +  'areapixapod' - Used in computation of sce fraction
%  +  'cone sce fraction' - SCE fraction for cone psfs
%
% Pupil and sointspread function
%  +  'pupil function' - The pupil function.  Must call wvfComputePupilFunction on wvf before get
%  +  'psf' - Point spread function.  Must call wvfComputePSF on wvf before get
%  +  'psf centered' - Peak of PSF is at center of returned matrix
%  +  '1d psf' - One dimensional horizontal (along row) slice through PSF centered on its max
%  +  'diffraction psf' - Diffraction limite PSF
%  +  'cone psf' - PSF as seen by cones for given weighting spectrum.
%
% Need to be implemented/checked/documented
%  +  'distanceperpix'
%  +  'samplesspace'
%  +  'strehl'     - Ratio of peak of diffraction limited to actual

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
%   7/20/12 dhb      Get rid of weighting spectrum, replace with cone psf info structure
%           dhb      Fix up areapix, areapixapod, scefrac and add comments.
%           dhb      Added more checking for stale pupil function and psf where needed.
%           dhb      Get of cone psf and cone sce fraction

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
        else
            widx = wvfWave2idx(wvf,varargin{1});
            val = wvf.zcoeffs(widx);
        end
        DIDAGET = true;
        
    case {'pupilsizemeasured','measuredpupilsize', 'measuredpupil', 'measuredpupilmm', 'measuredpupildiameter'}
        % Pupil diameter in mm over for which wavefront expansion is valid
        % wvfGet(wvf,'measured pupil','mm')
        % wvfGet(wvf,'measured pupil')
        val = wvf.measpupilMM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'wlmeasured','wavelengthmeasred','measuredwl', 'measuredwavelength'}
        % Measurement wavelength (nm)
        val = wvf.measWlNM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-9)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'opticalaxismeasued','measuredopticalaxis', 'measuredopticalaxisdeg'}
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
        
    case {'numberspatialsamples','spatialsamples', 'npixels', 'fieldsizepixels'}
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
        % wvfGet(wvf,'pupil plane size',units,wList)
        % Total size of computed field in pupil plane, for calculated
        % wavelengths(s)
                
        % Get wavelengths.  What if varargin{2} is empty?
        wList = varargin{2};
        waveIdx = wvfWave2idx(wvf,wList);
        wavelengths = wvfGet(wvf,'calc wavelengths','nm');
        
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
        
        % Unit conversion.  If varargin{1} is empty, then the units are
        % 'mm' and we leave it alone.
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        DIDAGET = true;
        
    case {'psfarcminpersample', 'psfarcminperpixel', 'arcminperpix'}
        % wvfGet(wvf,'psf arcmin per sample',wList)
        %
        % Arc minutes per pixel in psf domain, for calculated
        % wavelength(s).
        
        % Get wavelengths
        wavelengths = wvfGet(wvf,'calc wavelengths','mm');
        wList = varargin{1};
        waveIdx = wvfWave2idx(wvf,wList);
        
        % Figure out what's being held constant with wavelength and act
        % appropriately.
        whichDomain = wvfGet(wvf,'sample interval domain');
        if (strcmp(whichDomain,'psf'))
            val = wvfGet(wvf,'ref psf arcmin per pixel')*ones(length(waveIdx),1);
        elseif (strcmp(whichDomain,'pupil'))
            radiansPerPixel = ...
                wavelengths(waveIdx)/wvfGet(wvf,'ref pupil plane size','mm');
            val = (180*60/pi)*radiansPerPixel;
        else
            error('Unknown sample interval domain ''%s''',whichDomain);
        end
        DIDAGET = true;
        
    case {'psfanglepersample','angleperpixel','angperpix'}
        % Angular extent per pixel in the psf domain, for calculated
        % wavelength(s).
        %
        % wvfGet(wvf,'psf angle per sample',unit,wList)
        % unit = 'min' (default), 'deg', or 'sec'
        unit  = varargin{1};
        wList = varargin{2};
        val = wvfGet(wvf,'psf arcmin per sample',wList);
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
        
    case {'psfangularsamples','samplesangle','samplesarcmin','supportarcmin'}
        % Return one-d slice of sampled angles for psf, centered on 0, for
        % a single wavelength
        % wvfGet(wvf,'psf angular samples',unit,waveIdx)
        % unit = 'min' (default), 'deg', or 'sec'
        unit = varargin{1};
        wList = varargin{2};
        if (length(wList) > 1)
            error('This only works for one wavelength at a time');
        end
        anglePerPix = wvfGet(wvf,'psf angle per sample',unit,wList);
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'spatial samples');
        val = anglePerPix*((1:nPixels)-middleRow);
        DIDAGET = true;
        
    case {'psfspatialsamples','samplesspace','supportspace','spatialsupport'}
        % wvfGet(wvf,'samples space','um',wList)
        % Spatial support in samples, centered on 0
        % Unit and wavelength must be specified
        
        % This parameter matters for the OTF and PSF quite a bit.  It
        % is the number of um per degree on the retina.
        umPerDeg = (330*10^-6);
        
        unit = varargin{1}; wList = varargin{2};
        
        % Get the samples in degrees
        val = wvfGet(wvf,'psf angular samples','deg',wList);
        
        % Convert to meters and then to selected spatial scale
        val = val*umPerDeg;  % Sample in meters assuming 300 um / deg
        val = val*ieUnitScaleFactor(unit);
        DIDAGET = true;
        
    case {'pupilspatialsamples'}
        % wvfGet(wvf,'pupil spatial samples','mm',wList)
        % Spatial support in samples, centered on 0
        % Unit and wavelength must be specified
        
        unit = varargin{1}; wList = varargin{2};
        
        % Get the sampling rate in the pupil plane in space per sample
        spacePerSample = wvfGet(wvf,'pupil plane size',unit,wList)/wvfGet(wvf,'spatial samples');
        nSamples = wvfGet(wvf,'spatial samples');
        middleRow = wvfGet(wvf,'middle row');
        val = spacePerSample*((1:nSamples)-middleRow);
        DIDAGET = true;
       
    case {'middlerow'}
        val = floor(wvfGet(wvf,'npixels')/2) + 1;
        DIDAGET = true;
end

%% Calculation parameters
switch parm
    case {'calcpupilsize', 'calculatedpupil'}
        % Pupil size to assume when computing pupil function and PSF.  Must
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
        
    case {'calcwave','wave','calcwavelengths','wavelengths','wavelength','wls'}
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
        
    case {'calcconepsfinfo'}
        % Weighting spectrum used in calculation of polychromatic psf
        val = wvf.conePsfInfo;
        DIDAGET = true;
        
    case {'calcnwave','nwave','numbercalcwavelengths','nwavelengths'}
        % Number of wavelengths to calculate at
        val = length(wvf.wls);
        DIDAGET = true;
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
            wList = varargin{1};
            index = find(ismember(round(wave),round(wList)));
            if ~isempty(index), val = val(index);
            else error('Passed wavelength not contained in sceParams');
            end
        end
        DIDAGET = true;
        
    case {'scefraction','scefrac','stilescrawfordeffectfraction'}
        % How much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function before getting areapix.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        if isempty(varargin)
            vak = wvfGet(wvf,'area pixapod') ./ wvfGet(wvf,'areapix');
        else
            wList = varargin{1};
            val = wvfGet(wvf,'area pixapod',wList) ./ wvfGet(wvf,'areapix',wList);
        end
        DIDAGET = true;
        
    case {'areapix'}
        % This is the summed amplitude of the pupil function *before*
        % Stiles-Crawford correction over the pixels where the pupil
        % function is defined.  It doesn't have much physical significance,
        % but taking the ratio with areapixapod (just below) tells us
        % how much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function before getting areapix.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        if isempty(varargin)
            val = wvf.areapix;
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.areapix(idx);
            end
        end
        DIDAGET = true;
        
    case {'areapixapod'}
        % This is the summed amplitude of the pupil function *after*
        % Stiles-Crawford correction over the pixels where the pupil
        % function is defined.  It doesn't have much physical significance,
        % but taking the ratio with areapixapod (just above) tells us
        % how much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function before getting areapixapd.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        if isempty(varargin)
            val = wvf.areapixapod;
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.areapixapod(idx);
            end
        end
        DIDAGET = true;
        
        case {'conescefraction'}
        % SCE fraction for cone psfs

        % Can't do this unless psf is computed and not stale
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        [nil,coneSceFraction] = wvfComputeConePSF(wvf);
        DIDAGET = true;
end


% Computed pupil functions and point spread functions
switch parm
    case {'pupilfunction','pupilfunc','pupfun'}
        % The pupil function is derived from Zernicke coefficients in the
        % routine wvfComputePupilFunction If there are multiple
        % wavelengths, then this is a cell array of matrices
        %   wvfGet(wvf,'pupilfunc',wList)
        
        % Can't do the get unless it has already been computed and is not stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function before getting it.  Use wvfComputePupilFunction or wvfComputePSF.');
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
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.pupilfunc{idx};
            end
        end
        DIDAGET = true;
        
    case 'psf'
        % Get the ever-lovin' PSF.
        %   wvfGet(wvf,'psf',wList)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        % Return whole cell array of psfs over wavelength if
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
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.psf{idx};
            end
        end
        DIDAGET = true;
        
    case 'diffractionpsf'
        % Compute and return diffraction limited PSF.
        %   wvfGet(wvf,'diffraction psf',wList);
        if ~isempty(varargin), wList= varargin{1};
        else                   wList = wvfGet(wvf,'wave');
        end
        zcoeffs = zeros(65,1);
        wvfTemp = wvfSet(wvf,'zcoeffs',zcoeffs);
        wvfTemp = wvfSet(wvfTemp,'wave',wList(1));
        wvfTemp = wvfComputePSF(wvfTemp);
        val = wvfGet(wvfTemp,'psf',wList(1));
        DIDAGET = true;
        
    case 'strehl'
        % Strehl ratio. The strehl is the ratio of the peak of diff limited and the
        % existing PSF at each wavelength.
        %   wvfGet(wvf,'strehl',wList);
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        % We could write this so that with no arguments we return all of
        % the ratios across wavelengths.  For now, force a request for a
        % wavelength index.
        wList = varargin{1};
        psf = wvfGet(wvf,'psf',wList);
        dpsf = wvfGet(wvf,'diffraction psf',wList);
        val = max(psf(:))/max(dpsf(:));
        
        %         areaPixapod = wvfGet(wvf,'area pixapod',waveIdx);
        %         val = max(psf(:))/areaPixapod^2;
        %         % Old calculation was done in the compute pupil function routine.
        % Now, we do it on the fly in here, for a wavelength
        % strehl(wl) = max(max(psf{wl}))./(areapixapod(wl)^2);
        DIDAGET = true;
        
    case 'psfcentered'
        % PSF entered so that peak is at middle position in coordinate grid
        %   wvfGet(wvf,'psf centered',wList)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        if isempty(varargin), wList = wvfGet(wvf,'wave');
        else wList = varargin{1};
        end
        if length(wList) > 1, error('Only one wavelength permitted');
        else                  val = psfCenter(wvfGet(wvf,'psf',wList));
        end
        DIDAGET = true;
        
    case '1dpsf'
        % One dimensional slice through the PSF.
        %   wvfGet(wvf,'1d psf',wList,row)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        % Defaults
        wList = wvfGet(wvf,'wave');
        whichRow = wvfGet(wvf,'middle row');
        
        % Override with varargins
        if ~isempty(varargin),   wList    = varargin{1}; end
        if length(varargin) > 1, whichRow = varargin{2}; end
        
        psf = psfCenter(wvfGet(wvf,'psf',wList));
        val = psf(whichRow,:);
        DIDAGET = true;
        
    case 'conepsf'
        % PSF as seen by cones for specified weighting spectrum
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        % Defaults
        val = wvfComputeConePSF(wvf);

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

