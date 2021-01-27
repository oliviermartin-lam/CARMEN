close all; % close all figures
clear all; % clear the workspace
path_carmen = '/home/omartin/Projects/CARMEN/CARMEN/'; %% TBC %%
path_oomao = [path_carmen,'/_libOomao/']; % the oomao path
path_workspace = [path_carmen,'/_simu']; % the simulation folder path
addpath(genpath(path_oomao),genpath(path_carmen));

%%
flagMMSE = true; %% TBC %%
flagTraining = true; %% TBC %%
flagNoise = false; %% TBC %%

psTS = 0.219245;
S2Z =  (0.5 * 4.2 * 1e9 * psTS / (3600.*180./pi)) * fitsread('GLOB_mrz_7x7.fits')'; %in nm
    

%% DEFINING FIXED SIMULATION PARAMETERS
parFileCanary_3NGS
if ~flagNoise
    ron = 0;
end
%% DEFINING TOP-LEVEL CLASSES

%1\ ATMOSPHERE
atm   = atmosphere(photoNgs,r0,mean(L0),'layeredL0',L0,'fractionnalR0',fractionalR0,...
    'altitude',altitude,'windSpeed',windSpeed,'windDirection',windDirection);
%note that r0 is defined at phtoNgs wavelength

%2\ TELESCOPE
tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

%3\ TS SOURCE
sref = source('wavelength',photoNgs);

%4\ NGS
ngs      = source('magnitude',magNgs,'zenith',rNgs,'azimuth',dNgs,'wavelength',photoNgs);
nNgs    = numel(ngs);


%5\ TS WFS
% instantiating
ts = shackHartmann(nL,2*tel.resolution,minLightRatio);
ts.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
ts.camera.resolution        = [nPx*nL nPx*nL];

sref = sref.*tel*ts;
ts.INIT
+ts;
nSl = size(ts.slopes,1);
% Calibrate the pixel scale
ts.slopesUnits =  calibrateWfsPixelScale(ts,sref,wfsPscale*1e3,nPx);
tsPscale_simu      = ts.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx;
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',tsPscale_simu);
ts.camera.frameListener.Enabled = false;
ts.slopesListener.Enabled = false;
figure
imagesc(ts.camera)

%6\ NGS WFS
% instantiating
wfs = shackHartmann(nL,2*tel.resolution,minLightRatio);
wfs.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
wfs.camera.resolution        = [nPx*nL nPx*nL];
% Set up valid subaperture
sref = sref.*tel*wfs;
wfs.INIT
+wfs;

% Calibrate the pixel scale
wfs.slopesUnits =  calibrateWfsPixelScale(wfs,sref,wfsPscale*1e3,nPx);
wfsnPscale_simu      = wfs.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx;
wfs.camera.pixelScale = wfsnPscale_simu*constants.arcsec2radian;
fprintf('NGS WFS pixel scale set to %4.2f arcsec/pixel\n',wfsnPscale_simu);
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
wfs.camera.quantumEfficiency = wfsQE;

close all;

%% RUN THE END-TO-END SIMULATION
%Altitude layer height in km for the training set
nScreens = 15;
nIter = 1000;
alt = (1:1:nScreens)*1e3;
tmp1 = zeros(1,nScreens);
tmp2 = zeros(1,nScreens);
tmp3 = zeros(1,nScreens);
tmp4 = zeros(1,nScreens);

% Calculate the MMSE reconstructor optimized for a 7.5 km altitude layer
atm_mmse = atmosphere(photoNgs,r0,25,'fractionnalR0',[0.5,0.5],...
    'altitude',[0,2*1e3],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
Rmmse_2 = getMMSE_CANARY(tel,atm_mmse,ngs,wfs,sref,S2Z);
atm_mmse = atmosphere(photoNgs,r0,25,'fractionnalR0',[0.5,0.5],...
    'altitude',[0,7*1e3],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
Rmmse_7 = getMMSE_CANARY(tel,atm_mmse,ngs,wfs,sref,S2Z);
atm_mmse = atmosphere(photoNgs,r0,25,'fractionnalR0',[0.5,0.5],...
    'altitude',[0,12*1e3],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
Rmmse_12 = getMMSE_CANARY(tel,atm_mmse,ngs,wfs,sref,S2Z);

% loop on altitude
parfor l=1:15
    % GENERATE TELEMETRY
    atm = atmosphere(photoNgs,r0,25,'fractionnalR0',[0.5,0.5],...
        'altitude',[0,alt(l)],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
    trs = generateTelemetry(tel,atm,ngs,sref,wfs,ts,nIter,'training',flagTraining,'ron',0,'mmse',true,'S2Z',S2Z);
    
    % Get WAVEFRONT ERROR 
    
    tomoSl = Rmmse_2*reshape(trs.wfsSl,nSl*nNgs,nIter);
    tmp1(l) = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl - tomoSl);   
    
    tomoSl = Rmmse_7*reshape(trs.wfsSl,nSl*nNgs,nIter);
    tmp2(l) = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl - tomoSl);   
    
    tomoSl = Rmmse_12*reshape(trs.wfsSl,nSl*nNgs,nIter);
    tmp3(l) = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl - tomoSl);  
    
    tmp4(l) = trs.wfe.mmse;
    
end

wfe_mmse = [tmp1;tmp2;tmp3;tmp4];
fitswrite(wfe_mmse,'/home/omartin/Projects/CARMEN/RES/mmse_perf.fits')