function idx = wave2idx(wvf,wList)
% Convert the wavelength list (wList) to indices relative to
% wvfGet(wvf,'wave'); 
%
% For example, if wvfGet(wvf,'wave') is [400 500 600],
% and wList is [500,600] 
% then idx is [2,3].
%
% Should we return the nearest index or only exact matches?  We
% should start with only exact matches.
% Is this really only for 'wave'?  Should we have a flag for measured
% wavelength?
%
% Example
%  wvf = wvfCreate;
%  wvf = wvfSet(wvf,'wave',400:10:700);
%  wList = 500:100:700;
%  idx = wave2idx(wvf,wList)
%

wave = wvfGet(wvf,'wave');

% Check to within 1 nm
idx = find(ismember(round(wave),round(wList)));

if isempty(idx), warning('WVF:Wave2Idx','No matching wavelength %.02f',wList); end
return
