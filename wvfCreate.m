function wvfP = wvfCreate(varargin)
% Create the wavefront parameters structure
%
%    wvfP = wvfCreate(varargin)
%
% varargin:  Will be structured as param, val pairs
%   
% See also:  scePGet
%
% Examples:
%    wvfP = wvfCreate('wave',[400:10:700]);
%    wvfP = wvfCreate('
%
% (c) WVF Toolbox Team 2011

% Should have a switch for reading varargin param, val pairs.

% Wavelength samples - Default is 550, monochromatic
S = [550 1 1]; 

% Book-keeping
wvfP.name = 'default';
wvfP.type = 'wvf';

% Parameters in wvf
wvfP.zcoeffs     = zeros(65,1);     % Zernicke coefficients

% These are diameters. The measpupilMM is the pupil diameter when the
% aberrations were measured.  The reason for specifying the diameter at
% measurement time is to allow an error check that you're not trying to
% compute for a pupil larger than that for which you have data. (DHB)
wvfP.measpupilMM = 8;               % Maximum pupil diameter
wvfP.calcpupilMM = 3;               % Used for this calculation

%
wvfP.nominalFocusWl = 550;          % In focus wavelength (nm)
wvfP.defocusDiopters = 0;           % Defocus

% What field are we talking about here?
wvfP.sizeOfFieldPixels = 201;       % In pixels?  No units?
wvfP.sizeOfFieldMM = 16.212;        % Not sure which field

% Should this be just 1 wavelength or a vector? DHB seems to use wls as the
% wavelength list and S as a 3-vector of [start sample number]
wls = SToWls(S);
wvfP.wls         = wls;             % Wavelength samples

% Something about the cones.  
T = load('T_cones_ss2');   % Probably in the PTB
T_cones = SplineCmf(T.S_cones_ss2,T.T_cones_ss2,S);
% vcNewGraphWin; plot(wave,T_cones'); xlabel('Wavelength (nm)');

% Weighting spectrum, for combining the PSFs in an average 
T = load('spd_D65');
weightingSpectrum = SplineSpd(T.S_D65,T.spd_D65,S);

wvfP.T_cones = T_cones;             % Resampled cone spectral absorptions
wvfP.weightingSpectrum = weightingSpectrum;  % Probably used for combined psf

% Sets up the Stiles Crawford Effect parameters. 
wvfP.sceParams = sceCreate([],'none');

% Set additional arguments
if ~isempty(varargin)
    if isodd(length(varargin))
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        wvfP = wvfSet(wvfP,varargin{ii},varargin{ii+1});
    end
end

return
