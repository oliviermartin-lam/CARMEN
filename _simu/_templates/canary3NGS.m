%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
path_carmen = '/home/carlos/0_Projects/8_Conv_Carmen/matlab_sim/CARMEN/'; %% TBC %%
path_res = '/home/carlos/0_Projects/8_Conv_Carmen/matlab_sim/Data/New_code_test/'; %% TBC %%
path_oomao = [path_carmen,'/_libOomao/']; % the oomao path
path_workspace = [path_carmen,'/_simu']; % the simulation folder path
addpath(genpath(path_oomao),path_workspace);

flagMMSE = true; %% TBC %%
flagTraining = true; %% TBC %%
flagNoise = false; %% TBC %%
flagSave = true;

psTS = 0.219245;
S2Z =  (0.5 * 4.2 * 1e9 * psTS / (3600.*180./pi)) * fitsread('GLOB_mrz_7x7.fits')'; %in nm
    
% Creating saving folders
st = {'WFS_CAM','WFS_SLOPES','TS_CAM','TS_SLOPES','TOMO_SLOPES'};
for k=1:numel(st)
    if ~isfolder([path_res,st{k}])
        mkdir([path_res,st{k},'/'])
    end
end
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
%check where sources are in polar coordinates
figure;
pp = ngs.polar; 
for i=1:nNgs
    pp(i).MarkerFaceColor = 'b';
    pp(i).MarkerSize= 10;
    pp(i).Marker= 'o';
    pp(i).Color   = 'b';
end

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
if flagTraining
    alt = linspace(hmin/1e3,hmax/1e3,nScreens);    
end
%If flagTraining = false, the atmosphere is defined in parFileCanary_3NGS

for l=1:nScreens
    % GENERATE TELEMETRY
    if flagTraining
        %if training, update the altitude
        atm = atmosphere(photoNgs,r0,mean(L0_training),'layeredL0',L0_training,'fractionnalR0',fractionalR0_training,...
            'altitude',[0,alt(l)*1e3],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
    end
    trs = generateTelemetry(tel,atm,ngs,sref,wfs,ts,nIter,'training',flagTraining,'ron',ron,'mmse',flagMMSE,'S2Z',S2Z);

   % SAVE TELEMETRY    
   if flagSave
       fitswrite(trs.wfsSl,[path_res,'WFS_SLOPES/offaxiswfss_slopes_',num2str(atm.layer(2).altitude),'km.fits']);
       fitswrite(trs.tsSl,[path_res,'TS_SLOPES/ts_slopes_',num2str(atm.layer(2).altitude),'km.fits']);
       fitswrite(trs.wfsCam,[path_res,'WFS_CAM/offaxiswfss_cam_',num2str(atm.layer(2).altitude),'km.fits']);
       fitswrite(trs.tsCam,[path_res,'TS_CAM/ts_slopes_',num2str(atm.layer(2).altitude),'km.fits']);
       if flagMMSE
           fitswrite(trs.tomoSl,[path_res,'TOMO_SLOPES/tomo_slopes_',num2str(atm.layer(2).altitude),'km.fits']);
       end
   end
   trs.wfe
     
end