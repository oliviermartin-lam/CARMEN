close all; % close all figures
clear all; % clear the workspace

[~,msg] = unix('echo "$USER"');

if contains(msg,'omartin')
    path_carmen = '/home/omartin/Projects/CARMEN/CARMEN/'; %% TBC %%
    path_res    = '/home/omartin/Projects/CARMEN/RES/'; %% TBC %%
else
    path_carmen = '/home/carlos/0_Projects/8_Conv_Carmen/matlab_sim/CARMEN/'; %% TBC %%
    path_res    = '/home/carlos/0_Projects/8_Conv_Carmen/matlab_sim/Data/New_code_test/'; %% TBC %%
end

path_oomao      = [path_carmen,'_libOomao/']; % the oomao path
path_simu       = [path_carmen,'_simu']; % the simulation folder path
path_processing = [path_carmen,'_dataprocessing']; % the simulation folder path

addpath(genpath(path_oomao),genpath(path_simu),genpath(path_processing));

%% MANAGE OPTIONS
flagMMSE        = false; %% TBC %%
flagTraining    = true; %% TBC %%
flagNoise       = false; %% TBC %%
flagSave        = true;
simuCase        = 'Canary_3NGS'; % TBC


%% DEFINING FIXED SIMULATION PARAMETERS
% read the parameters file

if ~isfile([path_simu,'/_parFiles/','parFile',simuCase,'.m'])
   error('The parfile does not exist !');  
end
    
run(['parFile',simuCase])
if ~flagNoise
    ron = 0;
end
%1\ TELESCOPE
tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

%2\ TS SOURCE
sref = source('wavelength',photoGs);

%3\ GUIDE STARS
if hGs == 0
    % NATURAL GUIDE STAR
    gs      = source('magnitude',magGs,'zenith',rGs,'azimuth',dGs,'wavelength',photoGs);
else
    % LASER GUIDE STAR
    nPhoton = mean(photoGs.zeroPoint*10.^(-0.4*magGs));
    if nLayerNa > 0
        hNa     = linspace(hGs-naThickness/2,hGs+naThickness/2,nLayerNa);
    else
        naWeight= 1;
        hNa     = hGs;
    end
    gs      = laserGuideStar(tel.D/nL,tel.D,hGs,spotFwhm,nPhoton,naWeight,'zenith',rGs,'azimuth',dGs,'height',...
        hNa,'wavelength',photoGs,'viewPoint',viewPoint);
    
     % lgs = laserGuideStar(apertureSize,apertureDistance,...
    %        meanAltitude,fwhmInArcsec,nPhoton,naDensityProfile,sourceParams)
    %        creates a Laser Guide Star object with the size of the
    %        aperture used to image the LGS, the distance between the LGS
    %        launch telescope and the furthest aperture, the LGS mean
    %        altitude, the LGS intensity profil FHWM in arcsec, the number
    %        of photon per square meter and per second, the Sodium density
    %        altitude profile and the usual source parameters
    %
    % Example: A LGS on a 25m telescope with a 60x60 Shack-Hartmann WFS
    % launched from the telescope edge at an altitude of 90km and a fwhm of
    % 1 arcsec is created with:
    % lgs = laserGuideStar(25/60,25,90e3,1,1e6,ones(1,11),...
    %        'wavelength',photometry.Na,...
    %        'height',1e3*((-5:5)+90),...
    %        'viewPoint',[-25/2,0])
    % The number of photon is 1e6 and the Na density profile is flat.
    % The telescope pupil needs to be sampled with at least 60*44 pixels
    %
    % See also: source
    
end
nGs    = numel(gs);
%check where sources are in polar coordinates
figure;
pp = gs.polar; 
for i=1:nGs
    pp(i).MarkerFaceColor = 'b';
    pp(i).MarkerSize= 10;
    pp(i).Marker= 'o';
    pp(i).Color   = 'b';
end

                    
%4\ TS WFS
ts = shackHartmann(nL,tel.resolution,minLightRatio);
if ~nyquistFlag
    ts.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
    ts.camera.resolution        = [nPx*nL nPx*nL];
end

sref = sref.*tel*ts;
ts.INIT
+ts;
nSl = size(ts.slopes,1);
% Calibrate the pixel scale
ts.slopesUnits                  = calibrateWfsPixelScale(ts,sref,wfsPscale*1e3,nPx);
tsPscale_simu                   = ts.lenslets.fieldStopSize*constants.radian2arcsec*photoGs.wavelength/d/nPx;
ts.camera.frameListener.Enabled = false;
ts.slopesListener.Enabled       = false;
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',tsPscale_simu);

%5\ NGS WFS
% instantiating
wfs = shackHartmann(nL,tel.resolution,minLightRatio);
if ~nyquistFlag
    wfs.lenslets.fieldStopSize   = nPx*round(fovWfsinlod/nPx);
    wfs.camera.resolution        = [nPx*nL nPx*nL];    
end

% Set up valid subaperture
sref = sref.*tel*wfs;
wfs.INIT
+wfs;
% Calibrate the pixel scale
wfs.slopesUnits                     = calibrateWfsPixelScale(wfs,sref,wfsPscale*1e3,nPx);
wfsnPscale_simu                     = wfs.lenslets.fieldStopSize*constants.radian2arcsec*photoGs.wavelength/d/nPx;
wfs.camera.pixelScale               = wfsnPscale_simu*constants.arcsec2radian;
wfs.camera.frameListener.Enabled    = false;
wfs.slopesListener.Enabled          = false;
wfs.camera.quantumEfficiency        = wfsQE;
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',wfsnPscale_simu);

%6\ ATMOSPHERE
atm   = atmosphere(photoGs,r0,mean(L0),'layeredL0',L0,'fractionnalR0',fractionalR0,...
    'altitude',altitude,'windSpeed',windSpeed,'windDirection',windDirection);

%% ZERNIKE RECONSTRUCTION MATRIX
S2Z   = calibrateZernikeReconstructionMatrix(tel,ts,sref,'nModes',nZern,'amp',sref.wavelength/40);

%% GET PROFILES
config = 'variable';
nL_c = 4; % number of layers
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
for k=1:nDate
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
r0    = ((2*pi/photoGs.wavelength)^2 * 0.423*sum(Cn2_c,2)).^(-3/5);

figure;
plot(date_bin,r0);
ylabel('$r_0$ (m)','interpreter','latex','fontsize',18);
xlabel('Time (h)','interpreter','latex','fontsize',18);
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );

%% CREATING SAVING FOLDERS
if flagSave
    st = {'WFS_CAM','WFS_SLOPES','TS_CAM','TS_SLOPES','TOMO_SLOPES'};
    saveDir = [path_res,upper(simuCase),'/'];
    for k=1:numel(st)
        path = [saveDir,st{k}];
        if ~isfolder(path)
            mkdir(path)
        end
    end
end

%% GENERATE TELEMETRY

for kBin = 1:10
    
    %UPDATE ATMOSPHERE
    %note that we do not care about the tmeporal aspect as we're
    %going to generate independent phase screens
    
    atm   = atmosphere(photoGs,r0(kBin),25,'layeredL0',25,'fractionnalR0',Cn2_c(kBin,:)/sum(Cn2_c(kBin,:)),...
        'altitude',alt_c(kBin,:),'windSpeed',10*ones(1,nL_c),'windDirection',zeros(1,nL_c));
    
    % GENERATE TELEMETRY
    trs = generateTelemetry(tel,atm,gs,sref,wfs,ts,nIter,'training',true,'ron',ron,'mmse',false,'S2Z',S2Z);
    
    %MMSE RECONSTRUCTOR
    if kBin == 1
        Rmmse = getMMSE(tel,atm,gs,wfs,sref,S2Z);
    end
    trs.tomoSl = Rmmse*reshape(trs.wfsSl,nSl*nGs,nIter);
    
    % SAVE TELEMETRY
    fitswrite(trs.wfsSl,[saveDir,'WFS_SLOPES/''offaxiswfss_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tsSl,[saveDir,'TS_SLOPES/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.wfsCam,[saveDir,'WFS_CAM/offaxiswfss_cam_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tsCam,[saveDir,'TS_CAM/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
    fitswrite(trs.tomoSl,[saveDir,'TOMO_SLOPES/tomo_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',config,'_',num2str(date_bin(kBin)),'.fits']);
end