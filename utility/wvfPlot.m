function wvfPlot(wvfP,pType,varargin)
% Gateway routine for wavefront plots
%
%   wvfPlot(wvfP,pType,varargin);
%
% Angle units are 'sec','min', or 'deg' default - 'min'
% Space units are 'm','cm','mm','um' default - 'mm'
%
% Plot types:
%   1d psf angle - wvfPlot(wvfP,'1d psf angle',unit,maxVal);
%   2d psf angle - wvfPlot(wvfP,'1d psf angle',unit,maxVal);
%   1d psf space - wvfPlot(wvfP,'1d psf space',unit,maxVal);
%   2d psf space - wvfPlot(wvfP,'2d psf space',unit,maxVal);
%
% Examples
%
%
%
% (c) Wavefront Toolbox Team 2012

if ieNotDefined('wvfP'), error('Wavefront structure required.'); end
if ieNotDefined('pType'), pType = '1dpsf'; end

pType = ieParamFormat(pType);

switch(pType)
    case {'1dpsf','1dpsfangle','1dpsfnormalized'}
        % wvfPlot(wvfP,'1d psf angle',unit,plotRangeArcMin);
        % 
        unit = 'min';
        pRange = inf;  % Arc minutes
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, pRange = varargin{2}; end
        
        psfLine = wvfGet(wvfP,'1d psf');
        psfLine = psfLine/max(psfLine(:));
        samp = wvfGet(wvfP,'samples angle',unit);
        
        % Make a plot through of the returned PSF in the central region.
        index = find(abs(samp) < pRange);
        samp = samp(index);
        psfLine = psfLine(index);
        plot(samp,psfLine,'r','LineWidth',4);
        str = sprintf('Angle (%s)',unit);
        xlabel(str); ylabel('PSF slice')
        
        % Store the data
        udata.x = samp;
        udata.y = psfLine;
        set(gcf,'userdata',udata);
case {'2dpsf','2dpsfangle'}
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
        
        mesh(samp,samp,psf);
        xlabel('Angle'); ylabel('Angle'); zlabel('PSF')
        s = sprintf('Angle (%s)',unit); 
        xlabel(s); ylabel(s);
        zlabel('PSF amplitude')
        
        udata.x = samp;
        udata.y = samp;
        udata.z = psf;
        set(gcf,'userdata',udata);
    
case {'1dpsfspace','1dpsfspacenormalized'}
        % wvfPlot(wvfP,'1d psf normalized',plotRangeArcMin);
        % 
        unit = 'mm';
        pRange = inf;
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, pRange = varargin{2}; end
        
        psfLine = wvfGet(wvfP,'1d psf');
        psfLine = psfLine/max(psfLine(:));
        samp = wvfGet(wvfP,'spatial support',unit);
        
        % Make a plot through of the returned PSF in the central region.
        index = find(abs(samp) < pRange);
        samp = samp(index); psfLine = psfLine(index);
        plot(samp,psfLine,'r','LineWidth',4);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel('PSF slice')

        % Store the data
        udata.x = samp;
        udata.y = psfLine;
        set(gcf,'userdata',udata);

    case {'2dpsfspace'}
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
        
        % Extract within the range
        if length(varargin) > 1
         pRange = varargin{2}; 
         index = (abs(samp) < pRange);
         samp = samp(index);
         psf = psf(index,index);
        end
        
        mesh(samp,samp,psf);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel(s);
        zlabel('Relative amplitude')
        
        udata.x = samp;
        udata.y = samp;
        udata.z = psf;
        set(gcf,'userdata',udata);
        
    otherwise
        error('Unknown plot type %s\n',pType);
end

return
