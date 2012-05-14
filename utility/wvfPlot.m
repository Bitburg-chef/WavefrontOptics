function [uData, pData] = wvfPlot(wvfP,pType,varargin)
% Gateway routine for wavefront plots
%
%   [userData, plotData] = wvfPlot(wvfP,pType,varargin);
%
% Angle units are 'sec','min', or 'deg' default - 'min'
% Space units are 'm','cm','mm','um' default - 'mm'
%
% Plot types:
%  Mesh (default)
%   2d psf angle -   wvfPlot(wvfP,'2d psf angle',[unit], [waveIdx], [plotRangeArcMin]);
%   2d psf angle normalized -
%   2d psf space -   wvfPlot(wvfP,'2d psf space',[unit], [waveIdx], [plotRangeArcMin]);
%   2d psf space normalized -
%
%  Line
%   1d psf angle -
%   1d psf angle normalized - mesh
%   1d psf space -
%   1d psf space normalized -
%
%  PSF image
%   image psf - wvfPlot(wvfP,'image psf space',[unit],[waveIdx], [plotRange]);
%
%  Pupil FUnctions
%      wvfPlot(wvfP,'2d pupil amplitude space','mm',pRange)
%
%
%   2d pupil amplitude space (KP 3/11/12, in progress)
%   2d pupil phase space (MDL 3/18/12, in progress)
%   call with '2d pupil function space' for now. uses MM.
%
% Examples
%    [~,p]=wvfPlot(wvfP,'1d psf space',unit,maxVal);
%    set(p,'color','r')
%    [u,p]=wvfPlot(wvfP,'1d psf space',unit,maxVal);
%    figure; plot(u.x,u.y)
%
% Note: wvfPlot currently can only handle PSFs and Pupil Functions which
% are only calculated for 1 wavelength at a time (2d matrix instead of 3rd
% dimension of wavelengths).
% For wvfP with PSFs and Pupil Functions calculated over multiple
% wavelengths, it will be necessary to call wvfPlot on a separate wvf
% which only has PSF/pupilfunction for single wavelength.
% Ex: psf = wvfGet(wvfP,'psf'); %3d psf matrix with 3rd dim of wls
%     psf = psf(:,:,1) %pulls out 2d psf for 1 wavelength
%     wvftemp = wvfSet(wvfP,'psf',psf); %temp wvf with single wl psf
%     wvfPlot(wvftemp,'2dpsfspacenormalized','mm',3); %plots single wl
% KP 3/12/12
%
% (c) Wavefront Toolbox Team 2012

if ieNotDefined('wvfP'), error('Wavefront structure required.'); end
if ieNotDefined('pType'), pType = '1dpsf'; end

uData = [];
pType = ieParamFormat(pType);
theseArgs = [];

switch(pType)
    
    case {'2dpsf','2dpsfangle','2dpsfanglenormalized'}
        % wvfPlot(wvfP,'2d psf angle normalized',unit, waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples angle',unit);
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
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples space',unit);
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
        
        case {'imagepsfspace','imagepsfspacenormalized'}
        % wvfPlot(wvfP,'image psf space',unit,waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        samp = wvfGet(wvfP,'samples space',unit);
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
        
    case {'1dpsf','1dpsfangle','1dpsfanglenormalized'}
        % wvfPlot(wvfP,'1d psf angle',unit,waveIdx, plotRangeArcMin);
        %
        if ~isempty(varargin)
            [unit, waveIdx, pRange] = wvfReadArg(varargin);
        end
        
        psfLine = wvfGet(wvfP,'1d psf',waveIdx);
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
        
        samp = wvfGet(wvfP,'spatial support',unit);
        
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
        
    case {'2dpupilamplitudespace'}
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
        
        samp      = wvfGet(wvfP,'samples space');
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
        
    case {'2dpupilphasespace'}
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
        
        samp = wvfGet(wvfP,'samples space');
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



