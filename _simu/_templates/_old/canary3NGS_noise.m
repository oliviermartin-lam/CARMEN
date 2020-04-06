
%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
path_carmen = '/home/carlos/8_Conv_Carmen/matlab_sim/CODES'; %% TBC %%
path_res = '/home/carlos/8_Conv_Carmen/matlab_sim/Test_1/'; %% TBC %%
path_oomao = [path_carmen,'/_libOomao/']; % the oomao path
path_workspace = [path_carmen,'/_simu']; % the simulation folder path
addpath(genpath(path_oomao),path_workspace);

flagMMSE = true; %% TBC %%

% Creating saving folders
st = {'WFS_CAM','WFS_SLOPES','TS_CAM','TS_SLOPES','TOMO_SLOPES'};
for k=1:numel(st)
    if ~isfolder([path_res,st{k}])
        mkdir([path_res,st{k},'/'])
    end
end


%% DEFINING FIXED SIMULATION PARAMETERS
parFileCanary_3NGS

%% DEFINING TOP-LEVEL CLASSES

%1\ TELESCOPE
tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

%2\ TS SOURCE
sref = source('wavelength',photoNgs);

%3\ NGS
ngs      = source('magnitude',magNgs,'zenith',rNgs,'azimuth',dNgs,'wavelength',photoNgs);
nNgs    = numel(ngs);
%check where sources are in polar coordinates
figure;
pp = ngs.polar; 
for i=1:nNgs
    pp(i).MarkerFaceColor = 'b';
    pp(i).MarkerSize= 10;
    pp(i).Marker= 'o';
    pp(i).Color   = 'b';
end

%4\ TS WFS
% instantiating
ts = shackHartmann(nL,2*tel.resolution,minLightRatio);
ts.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
ts.camera.resolution        = [nPx*nL nPx*nL];

sref = sref.*tel*ts;
ts.INIT
+ts;
nSl = size(ts.slopes,1);
% Calibrate the pixel scale
ts.gainCalibration(tel,sref);
tsPscale_simu      = ts.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx;
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',tsPscale_simu);
ts.camera.frameListener.Enabled = false;
ts.slopesListener.Enabled = false;
figure
imagesc(ts.camera)

%5\ NGS WFS
% instantiating
wfs = shackHartmann(nL,2*tel.resolution,minLightRatio);
wfs.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
wfs.camera.resolution        = [nPx*nL nPx*nL];
% Set up valid subaperture
sref = sref.*tel*wfs;
wfs.INIT
+wfs;
% Calibrate the pixel scale
wfs.gainCalibration(tel,sref);
wfsnPscale_simu      = wfs.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx;
wfs.camera.pixelScale = wfsnPscale_simu*constants.arcsec2radian;
fprintf('NGS WFS pixel scale set to %4.2f arcsec/pixel\n',wfsnPscale_simu);
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

%% LOOP
dt    = 0;
count = 0;
nIter = round(nExp*Fe); %simulation time in iterations;
    
% Instantiate outputs
% Slopes time-vector
tsSl = zeros(nSl,nIter);
wfsSl = zeros(nSl,nNgs,nIter);
tomoSl = zeros(nSl,nIter);
% WFSs camera
tsCam = zeros(nPx*nL,nPx*nL,nIter);
wfsCam = zeros(nPx*nL,nL*nPx*nNgs,nIter);

%Enable readout and phton noise
wfs.camera.readOutNoise = ron;
wfs.camera.photonNoise = phNoise;
wfs.camera.quantumEfficiency = wfsQE;
    
% Loop on atmosphere status

alt = linspace(0.01,15.5,nScreens);

for l=1:nScreens
    
    % Redefining a telescope class to dissociate from the atmosphere
    tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

    % New atmosphere class
    altitude = [0, 4000., 10000., 15500.] %[0.,alt(l)]*1e3; %% TBC %%  
    atm   = atmosphere(photoAtm,r0,mean(L0),...
        'layeredL0',L0,'fractionnalR0',fractionalR0,...
        'altitude',altitude,'windSpeed',windSpeed,...
        'windDirection',windDirection);
    atm.wavelength = photoNgs;
    
    
    if flagMMSE       
        Rmmse = getMMSE_CANARY(tel,atm,ngs,wfs,sref);
    end
    
    
    %% RUN THE END-TO-END SIMULATION            
    tel = tel + atm; %Associate telescope and atmosphere;
    sref = sref.*tel*ts; %Define the light optical path for the Truth sensor;
    ngs = ngs.*tel*wfs; %Define the light optical path for the Truth sensor;
            
    close all;            
    for kIter=1:nIter
        tic;
        %1\ Updating phase screens
        +tel;
        
        %2\ Propagating sources to the TS
        +sref;
        %store TS slopes and pixels
        tsSl(:,kIter) = ts.slopes;
        tsCam(:,:,kIter) = ts.camera.frame;
        
        %3\ Propagating sources to the NGS WFS
        +ngs;
        %store slopes and pixels
        wfsSl(:,:,kIter) = wfs.slopes;
        wfsCam(:,:,kIter) = wfs.camera.frame;
        
        %4\ MMSE tomography
        if flagMMSE
            tomoSl(:,kIter) = Rmmse*wfs.slopes(:);
        end               
    end
    
    fitswrite(wfsSl,[path_res,'WFS_SLOPES/offaxiswfss_slopes_',num2str(alt(l)),'km.fits']);
    fitswrite(tsSl,[path_res,'TS_SLOPES/ts_slopes_',num2str(alt(l)),'km.fits']);
    fitswrite(wfsCam,[path_res,'WFS_CAM/offaxiswfss_cam_',num2str(alt(l)),'km.fits']);
    fitswrite(tsCam,[path_res,'TS_CAM/ts_slopes_',num2str(alt(l)),'km.fits']);
    if flagMMSE
    	fitswrite(tomoSl,[path_res,'TOMO_SLOPES/tomo_slopes_',num2str(alt(l)),'km.fits']);
    end
     
     getWaveFrontErrorFromSlopes_CANARY(tsSl)
end
