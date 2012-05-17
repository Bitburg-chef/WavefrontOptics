function [uData, pData] = wvfPlot(wvfP,pType,varargin)
% Wavefront plots
%
%   [userData, plotData] = wvfPlot(wvfP,pType,varargin);
%
% Plot types:
%   2d psf angle - mesh
%   2d psf space - mesh
%   2d OTF       - mesh (e.g., linepairs/'um')
%
%   1d psf angle - graph (middle horizontal line)
%   1d psf space - graph (middle horizontal line)
%
%   image psf    - image ('mm')
%   image pupil amp    - image ('mm')
%   image pupil phase  - image ('mm')
%
% Angle units are 'sec','min', or 'deg'   default - 'min'
% Space units are 'm','cm','mm','um'      default - 'mm'
%
% Add the string 'normalized' to force the 2d and 1d graphs to be scaled to
% a peak of 1.
%
% Examples
%    vcNewGraphWin; wvfP = wvfCreate; wvfP = wvfComputePSF(wvfP);
%    unit = 'um'; waveIdx = 1;
%    [u,p]= wvfPlot(wvfP,'1d psf space',unit,waveIdx);
%    set(p,'color','k','linewidth',2)
%
%    wvfPlot(wvfP,'image psf',unit,waveIdx);
%
% See also:  wvfComputePSF, wvfLoadHuman, vcNewGraphWin
%
% (c) Wavefront Toolbox Team 2012

if ieNotDefined('wvfP'), error('Wavefront structure required.'); end
if ieNotDefined('pType'), pType = '1dpsf'; end

uData = [];
pType = ieParamFormat(pType);

switch(pType)
    
    case {'2dpsf','2dpsfangle','2dpsfanglenormalized'}
        % wvfPlot(wvfP,'2d psf angle normalized',unit, waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples angle',unit,waveIdx);
        psf = wvfGet(wvfP,'psf',waveIdx);
        
        % Extract within the range
        if ~isempty(pRange)
            index = (abs(samp) < pRange);
            samp = samp(index);
            psf = psf(index,index);
        end
        
        % Search for key word normalized
        if ~isempty(strfind(pType,'normalized'))
            psf = psf(index,index)/max(psf(:));
        end
        
        % Start the plotting
        pData = mesh(samp,samp,psf);
        xlabel('Angle'); ylabel('Angle'); zlabel('PSF')
        s = sprintf('Angle (%s)',unit);
        xlabel(s); ylabel(s);
        zlabel('PSF amplitude')
        
        uData.x = samp; uData.y = samp; uData.z = psf;
        set(gcf,'userdata',uData);
        
    case {'2dpsfspace','2dpsfspacenormalized'}
        % wvfPlot(wvfP,'2d psf space',unit,waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples space',unit,waveIdx);
        psf = wvfGet(wvfP,'psf',waveIdx);
        if ~isempty(strfind(pType,'normalized'))
            psf = psf/max(psf(:));
        end
        
        % Extract within the range
        if ~isempty(pRange)
            index = (abs(samp) < pRange);
            samp = samp(index);
            psf = psf(index,index);
        end
                
        pData = mesh(samp,samp,psf);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel(s);
        zlabel('Relative amplitude')
        
        uData.x = samp; uData.y = samp; uData.z = psf;
        set(gcf,'userdata',uData);
        
        case {'imagepsf','imagepsfspace','imagepsfspacenormalized'}
        % wvfPlot(wvfP,'image psf space',unit,waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples space',unit,waveIdx);
        psf = wvfGet(wvfP,'psf',waveIdx);
        % If the string contains normalized
        if ~isempty(strfind(pType,'normalized'))
            psf = psf/max(psf(:));
        end
        
        % Extract within the range
        if ~isempty(pRange)
            index = (abs(samp) < pRange);
            samp = samp(index);
            psf = psf(index,index);
        end
          
        % Put up the image
        imagesc(samp,samp,psf); colormap(hot); axis image
        grid(gca,'on');
        set(gca,'xcolor',[.5 .5 .5]); set(gca,'ycolor',[.5 .5 .5]);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel(s);
        title('Relative amplitude')
        
        % Save the data
        uData.x = samp; uData.y = samp; uData.z = psf;
        set(gcf,'userdata',uData);
    
    case {'2dotf','otf'}
        % wvfPlot(wvfP,'2d otf',unit,waveIdx, plotRangeFreq);
        % wvfPlot(wvfP,'2d otf','mm',2, []);
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        % Get the data and if the string contains normalized ...
        psf = wvfGet(wvfP,'psf',waveIdx);
        if ~isempty(strfind(pType,'normalized'))
            psf = psf/max(psf(:));
        end
        
        % This stuff should move into wvfGet() - BW
        % Maybe freq = wvfGet(wvfP,'samples frequency',unit,waveIdx);
        samp = wvfGet(wvfP,'samples space',unit,waveIdx);
        nSamp = length(samp);
        dx = samp(2) - samp(1);
        nyquistF = 1 / (2*dx);   % Line pairs (cycles) per unit space
        freq = unitFrequencyList(nSamp)*nyquistF;

        % Compute OTF
        otf = fftshift(fft2(psf));
        
        % Restrict to parameter range
        if ~isempty(pRange)
            index = (abs(freq) < pRange);
            freq = freq(index);
            otf = otf(index,index);
        end
        
        % Axes, labeling, store data
        vcNewGraphWin; mesh(freq,freq,abs(otf))        
        str = sprintf('Freq (lines/%s)',unit);
        xlabel(str); ylabel(str);
        wave = wvfGet(wvfP,'wave','nm',waveIdx); 
        title(sprintf('OTF %.0f',wave));
        uData.fx = freq; uData.fy = freq; uData.otf = abs(otf);
        set(gcf,'userdata',uData);
        
    case {'1dpsf','1dpsfangle','1dpsfanglenormalized'}
        % wvfPlot(wvfP,'1d psf angle',unit,waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        psfLine = wvfGet(wvfP,'1d psf',waveIdx);
        samp = wvfGet(wvfP,'samples angle',unit,waveIdx);
        
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
        
    case {'1dpsfspace','1dpsfspacenormalized'}
        % wvfPlot(wvfP,'1d psf normalized',waveIdx,plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        psfLine = wvfGet(wvfP,'1d psf',waveIdx);
        if ~isempty(strfind(pType,'normalized'))
            psfLine = psfLine/max(psfLine(:));
        end
        
        samp = wvfGet(wvfP,'spatial support',unit, waveIdx);
        
        % Make a plot through of the returned PSF in the central region.
        if ~isempty(pRange)
            index = find(abs(samp) < pRange);
            samp = samp(index); psfLine = psfLine(index);
        end
        
        pData = plot(samp,psfLine,'r','LineWidth',4);
        s = sprintf('Position (%s)',unit);
        xlabel(s); ylabel('PSF slice')
        
        % Store the data
        uData.x = samp; uData.y = psfLine;
        set(gcf,'userdata',uData);
        
    case {'imagepupilamp','imagepupilampspace','2dpupilamplitudespace'}
        %wvfPlot(wvfP,'2d pupil amplitude space','mm',pRange)
        %plots the 2d pupil function amplitude for calculated pupil
        % Things to fix
        %  1. code in other plotting scales (distances or angles)
        %  2. fix units of pupil function plot
        
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        % need to change this, pupil function shouldn't be a mm related
        % plot...
        
        samp      = wvfGet(wvfP,'samples space','um',waveIdx);
        pupilfunc = wvfGet(wvfP,'pupil function',waveIdx);
        
        % Extract within the range
        if ~isempty(pRange)
            index = (abs(samp) < pRange);
            samp = samp(index);
            pupilfunc = pupilfunc(index,index);
        end
        
        pData = imagesc(samp,samp,abs(pupilfunc),[0 max(abs(pupilfunc(:)))]);
        s = sprintf('Position (%s)',unit);
        % this is a placeholder, need to fix with actual units?
        xlabel(s); ylabel(s);
        zlabel('Phase'); title('Pupil Function Amplitude'); colorbar;
        axis image;
        
    case {'imagepupilphase','2dpupilphasespace'}
        %plots the 2d pupil function PHASE for calculated pupil
        %
        %wvfPlot(wvfP,'2d pupil phase space','mm',pRange)
        %
        %some things to potentially fix:
        %1. modify colormap so that periodicity of phase is accounted for.
        %2. code in other plotting scales (distances or angles)
        %3. confirm plotting: currently 90deg flipped of wikipedia
        %4. somehow remove the 0 phase areas outside of calculated pupil
        %5. fix units of pupil function plot
        
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples space','um',waveIdx);
        pupilfunc = wvfGet(wvfP,'pupil function',waveIdx);
        
        % Extract within the range
        if ~isempty(pRange)
            index = (abs(samp) < pRange);
            samp = samp(index);
            pupilfunc = pupilfunc(index,index);
        end
        
        pData = imagesc(samp,samp,angle(pupilfunc),[-pi pi]);
        s = sprintf('Position (%s)',unit);
        % this is a placeholder, need to fix with actual units?
        xlabel(s); ylabel(s);
        zlabel('Phase'); title('Pupil Function Phase'); colorbar;
        axis image;
        
        
    otherwise
        error('Unknown plot type %s\n',pType);
end

return

end

%%% - Interpret the plotting arguments
function [units, waveIdx, pRange] = wvfReadArg(theseArgs)

if length(theseArgs) > 2, pRange = theseArgs{3};
else pRange = Inf;
end

if length(theseArgs) > 1, waveIdx = theseArgs{2};
else waveIdx = 1;
end

if ~isempty(theseArgs), units = theseArgs{1};
else units = 'min';
end

end



