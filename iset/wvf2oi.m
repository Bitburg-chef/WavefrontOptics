function oi = wvf2oi(wvfP, oType)
% Convert wavefront data to ISET optical image with optics
%
%  optics = wvf2oi(wvfP, [oType])
%
% Use Zernicke polynomial data in the wvfP structure and create a
% shift-invariant ISET optics model attached to the optical image
% structure.
%
% wvfP:  A wavefront parameters structure with a PSF
% oType: The type of oi structure 
%     Shift invariant -- Nothing special, default 'oi' and the wvfP data are
%     placed in the shift-invariant slot.  The optics type is set to shift
%     invariant 
%     
%     Human --  The optics is set up for the pupil size of the
%     wvfP structure, assuming a 17 mm focal length
%
%     Mouse --  Not yet implemented.
%
% Examples
%  pupilMM = 3; zCoefs = wvfLoadHuman(pupilMM);
%  wave = [400:10:700]';
%  wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
%  oi = wvf2oi(wvfP,'human');
%  fname = fullfile(isetRootPath,'data','optics','wvfHuman30.mat');
%  oi = oiSet(oi,'name','Human 3mm wvf');
%  vcExportObject(oi,fname);
%
% See also:
%
% Copyright Wavefront Toolbox Team 2012

%%
if ieNotDefined('oType'), oType = 'human'; end

% Create the shift-invariant PSF data structure
siData = wvf2PSF(wvfP);
pupil  = wvfGet(wvfP,'calculated pupil','m');

%% Create the OI
oType = ieParamFormat(oType);
switch oType
    case 'human'
        oi = oiCreate(oType);
        flength = 0.017;         % Human focal length is 17 mm
    
    case {'shiftinvariant','diffractionlimited','shift-invariant'}
        oi = oiCreate;
        flength = 0.017;         % Human focal length is 17 mm

    case 'mouse'
        flength = .003;          % Mouse focal length is 3 mm??
        error('Mouse not yet implemented');
        %         oi = oiCreate(oType);
    otherwise
        error('Unknown type %s\n',oType);
end

% Set up the optics and attach to OI
optics = siSynthetic('custom',oi,siData);
optics = opticsSet(optics,'model','shiftInvariant');
optics = opticsSet(optics,'fnumber',flength/pupil);
optics = opticsSet(optics,'flength',flength);
oi     = oiSet(oi,'optics',optics);

return
