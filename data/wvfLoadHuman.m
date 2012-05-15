function [sample_mean sample_cov] = vwfLoadHuman(pupilDiameterMM) %#ok<STOUT>
% Load  mean and covariance of Zernicke coefficients for human data
%
%   [sample_mean sample_cov] = vwfLoadHuman(pupilDiameterMM)
%
% The data are from Autrusseau et al. XXX
% Say more more.
%
% Example:
%   pd = 6.0;
%   [sample_mean sample_cov] = vwfLoadHuman(pd);
%
% See also: s_ThiboxModel
%
% Copyright Wavefront Toolbox Team, 2012

if ieNotDefined('pupilDiameterMM'), pupilDiameterMM = 6; end

%%
switch pupilDiameterMM
    case 7.5
        load('IASstats75','S','sample_mean');
    case 6.0
        load('IASstats60','S','sample_mean');
    case 4.5
        load('IASstats45','S','sample_mean');
    case 3.0
        load('IASstats30','S','sample_mean');
        
    otherwise
        error('Unknown pupil size')
end

sample_cov = S;

end