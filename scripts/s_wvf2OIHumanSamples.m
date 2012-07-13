%% s_wvf2OIHuman
%
% Check the variation in the Thibos wavefront model.
%
% See also: s_wvf2OIHuman and s_wvf2OI
%
% (c) Wavefront Toolbox Team, 2012


%%
s_initISET
maxUM = 10;
wave = 400:10:700; wave = wave(:);
pupilMM = 3; 

% Create some examples
[sample_mean S] = wvfLoadHuman(pupilMM);
N = 10;
zSamples = ieMvnrnd(sample_mean,S,N)';  
nCoeffs = size(zSamples,1);

%% Convert WVF human data to ISET
oiD = cell(N,1);

% Create samples
for ii=1:N
    name = sprintf('%d human-%d',ii,pupilMM);
    wvfP = wvfCreate('wave',wave,'name',name);
    z = wvfGet(wvfP,'zcoeffs');
    z(1:nCoeffs) = zSamples(:,ii);
    wvfP = wvfComputePSF(wvfP);
    oiD{ii} = wvf2oi(wvfP,'human');
    oiD{ii} = oiSet(oiD{ii},'name',name);
end

%% Now compare the slanted bar response in the OI 
% These are reasonably close for calculations separated by so many years.

scene  = sceneCreate('slanted bar');
scene  = sceneSet(scene,'h fov',1);
bb = blackbody(sceneGet(scene,'wave'),6500,'energy');
scene = sceneAdjustIlluminant(scene,bb);
% vcAddAndSelectObject(scene); sceneWindow;

sensor = sensorCreate('human');
sensor = sensorSet(sensor,'exp time',0.050);
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'hfov'),scene,oiD{1});

%%
uData = cell(N,1);
for ii=1:N
    oiD{ii} = oiCompute(oiD{ii},scene);
    sensorD = sensorCompute(sensor,oiD{ii});
    sensorD = sensorSet(sensorD,'name','Thibos calc');
    vcAddAndSelectObject(sensorD); sensorWindow('scale',1);
    [uData{ii},g] = plotSensor(sensorD,'electrons hline',[1,80]);
    close(g)
end

%%
vcNewGraphWin([],'tall');
hold on
c = {'r-','g-','b-'};
for ii=1:N
    for jj=1:3
        hold on; grid on;
        subplot(3,1,jj)
        plot(uData{ii}.pos{jj},uData{ii}.data{jj},c{jj})     
    end
end

%%
% sensorMW = sensorCompute(sensor,oiMW);
% sensorMW = sensorSet(sensorMW,'name','MW calc');
% vcAddAndSelectObject(sensorMW); sensorWindow('scale',1);
% plotSensor(sensorMW,'electrons hline',[1,80]);
% title('Marimont and Wandell')

%% End
