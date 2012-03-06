function [uData, pData] = wvfPlot(wvfP,pType,varargin)
% Gateway routine for wavefront plots
%
%   [userData, plotData] = wvfPlot(wvfP,pType,varargin);
%
% Angle units are 'sec','min', or 'deg' default - 'min'
% Space units are 'm','cm','mm','um' default - 'mm'
%
% Plot types:
%   1d psf angle - 
%   1d psf angle normalized - 
%   2d psf angle - 
%   2d psf angle normalized -
%   1d psf space - 
%   1d psf space normalized - 
%   2d psf space - 
%   2d psf space normalized - 
%
% Examples
%    [~,p]=wvfPlot(wvfP,'1d psf space',unit,maxVal);
%    set(p,'color','r')
%    [u,p]=wvfPlot(wvfP,'1d psf space',unit,maxVal);
%    figure; plot(u.x,u.y)
%
% (c) Wavefront Toolbox Team 2012

if ieNotDefined('wvfP'), error('Wavefront structure required.'); end
if ieNotDefined('pType'), pType = '1dpsf'; end

uData = [];
pType = ieParamFormat(pType);

switch(pType)
    case {'1dpsf','1dpsfangle','1dpsfanglenormalized'}
        % wvfPlot(wvfP,'1d psf angle',unit,plotRangeArcMin);
        % 
        unit = 'min';
        pRange = inf;  % Arc minutes
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, pRange = varargin{2}; end
        
        psfLine = wvfGet(wvfP,'1d psf');
        % psfLine = psfLine/max(psfLine(:));
        samp = wvfGet(wvfP,'samples angle',unit);
        
        % Make a plot through of the returned PSF in the central region.
        index = find(abs(samp) < pRange);
        samp = samp(index);
        psfLine = psfLine(index);
        if ~isempty(strfind(pType,'normalized'))
            psfLine = psfLine/max(psfLine(:));
        end
        
        pData = plot(samp,psfLine,'r','LineWidth',4);
        str = sprintf('Angle (%s)',unit);
        xlabel(str); ylabel('PSF slice')
        
        % Store the data
        uData.x = samp; uData.y = psfLine;
        set(gcf,'userdata',uData);
        
case {'2dpsf','2dpsfangle','2dpsfanglenormalized'}
        % Plot of the psf as a mesh, within a range
        % Samples are in arc minutes
        %   wvfPlot(wvfP,'2d psf',unit,pRange)
        %   wvfPlot(wvfP,'2d psf','min',2);
        %
        if isempty(varargin), unit = 'min'; 
        else unit = varargin{1};
        end
        samp = wvfGet(wvfP,'samples angle',unit);
        psf = wvfGet(wvfP,'psf');

        
        % Extract within the range
        if length(varargin) > 1
         pRange = varargin{2}; 
         index = (abs(samp) < pRange);
         samp = samp(index);
         psf = psf(index,index);
        end
        
        if ~isempty(strfind(pType,'normalized'))
            psf = psf(index,index)/max(psf(:));
        end
        pData = mesh(samp,samp,psf);
        xlabel('Angle'); ylabel('Angle'); zlabel('PSF')
        s = sprintf('Angle (%s)',unit); 
        xlabel(s); ylabel(s);
        zlabel('PSF amplitude')
        
        uData.x = samp; uData.y = samp; uData.z = psf;
        set(gcf,'userdata',uData);
    
case {'1dpsfspace','1dpsfspacenormalized'}
        % wvfPlot(wvfP,'1d psf normalized',plotRangeArcMin);
        % 
        unit = 'mm';
        pRange = inf;
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, pRange = varargin{2}; end
        
        psfLine = wvfGet(wvfP,'1d psf');
        if ~isempty(strfind(pType,'normalized'))
            psfLine = psfLine/max(psfLine(:)); 
        end

        % 
        samp = wvfGet(wvfP,'spatial support',unit);
        
        % Make a plot through of the returned PSF in the central region.
        index = find(abs(samp) < pRange);
        samp = samp(index); psfLine = psfLine(index);
        pData = plot(samp,psfLine,'r','LineWidth',4);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel('PSF slice')

        % Store the data
        uData.x = samp; uData.y = psfLine;
        set(gcf,'userdata',uData);

    case {'2dpsfspace','2dpsfspacenormalized'}
        % Plot of the psf as a mesh, within a range
        %
        % wvfPlot(wvfP,'2d psf',unit,pRange)
        % wvfPlot(wvfP,'2d psf','mm',2);
        %
        if isempty(varargin), unit = 'mm';
        else unit = varargin{1};
        end
        
        samp = wvfGet(wvfP,'samples space',unit);
        psf = wvfGet(wvfP,'psf');
        if ~isempty(strfind(pType,'normalized'))
            psf = psf/max(psf(:)); 
        end

        % Extract within the range
        if length(varargin) > 1
         pRange = varargin{2}; 
         index = (abs(samp) < pRange);
         samp = samp(index);
         psf = psf(index,index);
        end
        
        % Scale to one or not
        
        pData = mesh(samp,samp,psf);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel(s);
        zlabel('Relative amplitude')
        
        uData.x = samp; uData.y = samp; uData.z = psf;
        set(gcf,'userdata',uData);

    otherwise
        error('Unknown plot type %s\n',pType);
end

return
