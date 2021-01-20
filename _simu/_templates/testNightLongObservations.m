close all;clear all;clc;
path_carmen     = '/home/omartin/Projects/CARMEN/CARMEN/'; %% TBC %%
path_oomao      = [path_carmen,'/_libOomao/']; % the oomao path
path_workspace  = [path_carmen,'/_simu']; % the simulation folder path
addpath(genpath(path_oomao),genpath(path_carmen));

%% SAVING FOLDERS
path_res = '/home/omartin/Projects/CARMEN/RES/TIMESERIES/'; %TBC

% Creating saving folders
st = {'WFS_CAM','WFS_SLOPES','TS_CAM','TS_SLOPES','TOMO_SLOPES'};
for k=1:numel(st)
    if ~isfolder([path_res,st{k}])
        mkdir([path_res,st{k},'/'])
    end
end

%% DEFINING FIXED SIMULATION PARAMETERS
parFileCanary_3NGS

psTS = 0.219245;
S2Z =  (0.5 * 4.2 * 1e9 * psTS / (3600.*180./pi)) * fitsread('GLOB_mrz_7x7.fits')'; %in nm

%1\ TELESCOPE
tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

%2\ TS SOURCE
sref = source('wavelength',photoNgs);

%3\ NGS
ngs      = source('magnitude',magNgs,'zenith',rNgs,'azimuth',dNgs,'wavelength',photoNgs);
nNgs    = numel(ngs);

%4\ TS WFS
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
wfs.slopesUnits =  calibrateWfsPixelScale(wfs,sref,wfsPscale*1e3,nPx);
wfsnPscale_simu      = wfs.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx;
wfs.camera.pixelScale = wfsnPscale_simu*constants.arcsec2radian;
fprintf('NGS WFS pixel scale set to %4.2f arcsec/pixel\n',wfsnPscale_simu);
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
wfs.camera.quantumEfficiency = wfsQE;
close all;

%% GET PROFILES
config = 'typical';
nL_c = 7; % number of layers
dt   = 5; %in minutes

path_profiles = '/home/omartin/Projects/CARMEN/DOCUMENTS/profiles/';

switch config
    case 'calm'
        file = '20140716.txt';
    case 'typical'
        file = '20141006.txt';
    case 'variable'
        file = '20151005.txt';
end

%For each night file there is a row per profile with format:
%YYYY-MM-DDTHH:MM:SS.ss,r0,seeing,coherenceTime,isoplanaticAngle,scintillationIndex,alt_0,cn2_0,windSpeed_0,windDirection_1,al_t1,cn2_1, windSpeed _1, windDirection _1….alt_n,cn2_n, windSpeed _n, windDirection _n

% read
T     = readtable([path_profiles,file],'HeaderLines',0,'Delimiter',',');
nDate = numel(T{:,1});
alt   = 0:250:24750;
nL    = numel(alt);
Cn2   = zeros(nDate,nL);
loc_date = zeros(1,nDate);
for k=1:301
    tmp = cell2mat(T{k,1});
    loc_date(k) = str2double(tmp(12:13)) + str2double(tmp(15:16))/60 + str2double(tmp(18:19))/3600;
    if loc_date(k) > 12
        loc_date(k) = loc_date(k)-24;
    end
    for j=1:nL
        Cn2(k,j)  = T{k,8 + 4*(j-1)};
    end
end


% Temporal decimation
dt_ori = abs(median(diff(loc_date)))*60;
n_dec  = max(1,round(dt/dt_ori));
Cn2_bin = Cn2(1:n_dec:end,:);
date_bin = loc_date(1:n_dec:end);
% Compress
nBin = size(Cn2_bin,1);
Cn2_c = zeros(nBin,nL_c);
alt_c = zeros(nBin,nL_c);
for k=1:nBin
    [Cn2_c(k,:),alt_c(k,:)] = eqLayers(Cn2_bin(k,:),alt,nL_c);
end
r0    = ((2*pi/photoNgs.wavelength)^2 * 0.423*sum(Cn2_c,2)).^(-3/5);

figure;
plot(date_bin,r0);
ylabel('$r_0$ (m)','interpreter','latex','fontsize',18);
xlabel('Time (h)','interpreter','latex','fontsize',18);
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );

%% GENERATE TELEMETRY

for kBin = 1:1
    
    %UPDATE ATMOSPHERE
    %note that we do not care about the tmeporal aspect as we're
    %going to generate independent phase screens
    
    atm   = atmosphere(photoNgs,r0(kBin),25,'layeredL0',25,'fractionnalR0',Cn2_c(kBin,:)/sum(Cn2_c(kBin,:)),...
        'altitude',alt_c(kBin,:),'windSpeed',10*ones(1,nL_c),'windDirection',zeros(1,nL_c));
    
    % GENERATE TELEMETRY
    trs = generateTelemetry(tel,atm,ngs,sref,wfs,ts,nIter,'training',true,'ron',ron,'mmse',false,'S2Z',S2Z);
    
    %MMSE RECONSTRUCTOR
    if kBin == 1
        Rmmse = getMMSE_CANARY(tel,atm,ngs,wfs,sref,S2Z);
    end
    trs.tomoSl = Rmmse*reshape(trs.wfsSl,nSl*nNgs,nIter);
    
    % SAVE TELEMETRY
    fitswrite(trs.wfsSl,[path_res,'WFS_SLOPES/''offaxiswfss_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tsSl,[path_res,'TS_SLOPES/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.wfsCam,[path_res,'WFS_CAM/offaxiswfss_cam_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tsCam,[path_res,'TS_CAM/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tomoSl,[path_res,'TOMO_SLOPES/tomo_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
end