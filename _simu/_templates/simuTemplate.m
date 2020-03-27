
%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
path_carmen = '/home/omartin/Projects/CARMEN/CODES';
path_oomao = [path_carmen,'/_libOomao/']; % the oomao path
path_workspace = [path_carmen,'/_simu']; % the simulation folder path
addpath(genpath(path_oomao),path_workspace);

%% DEFINING SIMULATION PARAMETERS

%1\ TELESCOPE
D = 4.2; %telescope diameter in meter
cobs = 0.285; % central obscuration ratio
Fe       = 150; % Phase screen upgrate rate in Hz

%2\ ATMOSPHERE
photoAtm         = photometry.V0;      % photometric band to define atmospheric parameters V0->500nm
r0                     = 0.15;                         % coherence lenght in meters at 500 nm
fractionalR0      = [65,10,15,10]/100.; % fractional weight of turbulent layers
altitude             = [0.,4.,10.,16.]*1e3;  % altitude vector in meters
windSpeed        = [8.,12.,15.,18.];       % wind speed vector in m/s
windDirection    = [0,pi/2,pi/3,pi/4];    % wind direction vector in rad
L0                     = [25 25 25 25];         % Outer scale in meters

%3\ LGSs
photoLgs          = photometry.Rayleigh;                % photometric band    532 nm
xLgs                 = [-16.2635,16.2635,-16.2635,16.2635];                    % LGSs x cartesian position in arcsec
yLgs                 = [16.2635,16.2635,-16.2635,-16.2635];                    % LGSs y cartesian position in arcsec
[zS,aS]     = utilities.arcsec2polar(xLgs,yLgs);
hLgs                 = 21e3; %LGSs height
magLgs             = 8;

%4\ LGS WFS
% wfs geometry
nL                   = 7; %\# 1D lenselts
nPx                 = 8; %1D \# pixels per lenslet
nPxWfs           = nL*nPx;
d                     = D/nL; %subaperture size
minLightRatio  = 0.5; % ratio of subap illumination to be valid
wfsPscale         = 0.47;% Pixel scale in arcsec
fovWfs             = wfsPscale*nPx; %WFS field of view in arcsec
% Calcluate the number of pixels to get the closest pixel scale
lambdaOverd = constants.radian2arcsec*photoLgs.wavelength/d;
fovWfsinlod   = floor(fovWfs/lambdaOverd); % WFS fov in l/d units
Samp            = lambdaOverd/wfsPscale; % WFS sampling;
resTel           = nPx*nL*round(1/Samp); % \# pixels within the pupil
% wfs detector
brightestPixel = 20; %\# of brightest pixels used for slopes computation
ron = 4; %ron in e-
QE = 0.7; % Detector quantum efficiency
T = 0.1; % system transmittance
skyBg = 10; % number of background photons
darkBg = 0.2; % dark background
pixelGains = 1e3; % Gain

%6\ NGSs
photoNgs          = photometry.R;                % photometric band    640 nm
rNgs = [40.6 53 47.9]*constants.arcsec2radian;
dNgs = [pi/4 pi/2+pi/6,pi+pi/3];
magNgs             = [10.2 8.7 9.9];

%7\ NGS WFS
% wfs geometry
nL_n                = 7; %\# 1D lenselts
d_n                  = D/nL_n; %subaperture size
nPx_n              = 16; %1D \# pixels per lenslet
wfsPscale_n     = 0.26;% Pixel scale in arcsec
fovWfs_n         = wfsPscale_n*nPx_n; %WFS field of view in arcsec
% Calcluate the number of pixels to get the closest pixel scale
lambdaOverd_n = constants.radian2arcsec*photoNgs.wavelength/d_n;
fovWfsninlod   = floor(fovWfs_n/lambdaOverd_n); % WFS fov in l/d units
Samp_n            = lambdaOverd_n/wfsPscale_n; % WFS sampling;
% wfs detector
brightestPixel_n = 20; %\# of brightest pixels used for slopes computation
ron_n = 4; %ron in e-
QE_n = 0.6; % Detector quantum efficiency
T_n = 0.1; % system transmittance
skyBg_n = 100; % number of background photons
darkBg_n = 0.4; % dark background
pixelGains_n = 1e3; % Gain

fovTel = 2*max(rNgs)*constants.radian2arcsec; % telescope fov ni arcsec.

%% DEFINING TOP-LEVEL CLASSES

%1\ TELESCOPE
tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs,...
    'fieldOfViewInArcsec',fovTel,'samplingTime',1/Fe);

%2\ ATMOSPHERE
atm   = atmosphere(photoAtm,r0,mean(L0),...
    'layeredL0',L0,'fractionnalR0',fractionalR0,...
    'altitude',altitude,'windSpeed',windSpeed,...
    'windDirection',windDirection);
atm.wavelength = photoLgs;

%3\ TS SOURCE
sref = source('wavelength',photoLgs);

%4\ LGS
lgs      = source('magnitude',magLgs,'zenith',zS,'azimuth',aS,'wavelength',photoLgs,'height',hLgs);
nLgs    = numel(lgs);

%5\ NGS
ngs      = source('magnitude',magNgs,'zenith',rNgs,'azimuth',dNgs,'wavelength',photoNgs);
nNgs    = numel(ngs);
%check where sources are in polar coordinates
gs = [ngs lgs];
figure;
pp = gs.polar; 
for i=1:nNgs
    pp(i).MarkerFaceColor = 'b';
    pp(i).MarkerSize= 10;
    pp(i).Marker= 'o';
    pp(i).Color   = 'b';
end
for i=nNgs+1:nNgs +nLgs
    pp(i).MarkerFaceColor = 'g';
    pp(i).MarkerSize= 10;
    pp(i).Marker= 'd';
    pp(i).Color= 'g';
end

%5\ TS WFS
% instantiating
ts = shackHartmann(nL,2*tel.resolution,minLightRatio);
ts.lenslets.fieldStopSize   = 20;
ts.camera.resolution        = [nPx*nL nPx*nL];
% Set up valid subaperture
sref = sref.*tel*ts;
ts.INIT
+ts;
nSl = size(ts.slopes,1);
% Calibrate the pixel scale
ts.gainCalibration(tel,sref);
tsPscale_simu      = ts.lenslets.fieldStopSize*constants.radian2arcsec*photoLgs.wavelength/d/nPx;
%ts.slopesUnits = calibrateWfsPixelScale(ts,sref,tsPscale_simu*1e3,nPx);
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',tsPscale_simu);
ts.camera.frameListener.Enabled = false;
ts.slopesListener.Enabled = false;
figure
imagesc(ts.camera)

%6\ LGS WFS
% instantiating
wfs = shackHartmann(nL,2*tel.resolution,minLightRatio);
wfs.lenslets.fieldStopSize   = fovWfsinlod;
wfs.camera.resolution        = [nPx*nL nPx*nL];
% Set up valid subaperture
sref = sref.*tel*wfs;
wfs.INIT
+wfs;
% Calibrate the pixel scale
wfs.gainCalibration(tel,sref);
wfsPscale_simu      = wfs.lenslets.fieldStopSize*constants.radian2arcsec*photoLgs.wavelength/d/nPx;
wfs.camera.pixelScale = wfsPscale_simu*constants.arcsec2radian;
%wfs.slopesUnits = calibrateWfsPixelScale(wfs,sref,wfsPscale_simu*1e3,nPx);
fprintf('WFS pixel scale set to %4.2f arcsec/pixel\n',wfsPscale_simu);
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
% WFS noise
wfs.camera.photonNoise = true; %photon noise disactivated
wfs.camera.readOutNoise = ron; %ron in e-
wfs.camera.quantumEfficiency = QE*T;
wfs.camera.nPhotonBackground = skyBg;
wfs.camera.darkBackground = darkBg;
wfs.camera.pixelGains = pixelGains;
wfs.brightestPixel = brightestPixel;
% Noise variance; to be checked !
varn_lgs = mean(mean(wfs.theoreticalNoise(tel,atm,lgs,sref,'skyBackground',skyBg,'emccd',1,'centroidingAlgorithm','cog'))); %in phase rd^2
varn_lgs = varn_lgs*(constants.radian2arcsec*photoLgs.wavelength/2/pi/d/Samp^2)^2;


%7\ NGS WFS
% instantiating
wfsn = shackHartmann(nL_n,2*tel.resolution,minLightRatio);
wfsn.lenslets.fieldStopSize   = nPx_n*round(fovWfsninlod/nPx_n);
wfsn.camera.resolution        = [nPx_n*nL_n nPx_n*nL_n];
% Set up valid subaperture
sref = sref.*tel*wfsn;
wfsn.INIT
+wfsn;
% Calibrate the pixel scale
wfsn.gainCalibration(tel,sref);
wfsnPscale_simu      = wfsn.lenslets.fieldStopSize*constants.radian2arcsec*photoNgs.wavelength/d/nPx_n;
wfsn.camera.pixelScale = wfsnPscale_simu*constants.arcsec2radian;
%wfsn.slopesUnits = calibrateWfsPixelScale(wfsn,sref,wfsnPscale_simu*1e3,nPx_n);
fprintf('NGS WFS pixel scale set to %4.2f arcsec/pixel\n',wfsnPscale_simu);
wfsn.camera.frameListener.Enabled = false;
wfsn.slopesListener.Enabled = false;
% NGS WFS noise
wfsn.camera.photonNoise = true; %photon noise disactivated
wfsn.camera.readOutNoise = ron_n; %ron in e-
wfsn.camera.quantumEfficiency = QE_n*T_n;
wfsn.camera.nPhotonBackground = skyBg_n;
wfsn.camera.darkBackground = darkBg_n;
wfsn.camera.pixelGains = pixelGains_n;
wfsn.brightestPixel = brightestPixel_n;
% Noise variance; to be checked !
varn_ngs = mean(mean(wfsn.theoreticalNoise(tel,atm,ngs,sref,'skyBackground',skyBg_n,'emccd',1,'centroidingAlgorithm','cog'))); %in phase rd^2
varn_ngs = varn_ngs*(constants.radian2arcsec*photoNgs.wavelength/2/pi/d/Samp_n^2)^2;


%8\ JOINT MMSE TOMOGRAPHIC RECONSTRUCTOR
covModel = slopesCovarianceModel(tel,atm,wfs,[sref,ngs,lgs],'wfsType',[1,3,3,3,2,2,2,2],'noiseVar',[0,0,varn_ngs*ones(1,2*nNgs),varn_lgs*ones(1,2*nLgs)]);
covMat = covModel.getCovarianceMatrix();
R_LA = covModel.getMMMSEreconstructor(covMat); % MMSE reconstructor NGS+LGS slopes -> TS slopes

%9\ LGS MMSE TOMOGRAPHIC RECONSTRUCTOR
covModel_LGS = slopesCovarianceModel(tel,atm,wfs,[sref,lgs],'wfsType',[1,2,2,2,2],'noiseVar',[0,0,varn_lgs*ones(1,2*nLgs)]);
covMat_LGS = covModel_LGS.getCovarianceMatrix();
R_LGS = covModel_LGS.getMMMSEreconstructor(covMat_LGS); % MMSE reconstructor LGS slopes -> tip-tilt excluded TS slopes
% Tip-tilt matrix
t  = ones(wfsn.nValidLenslet,1);
tt = [t 0*t; 0*t,t];
R_TT = pinv(eye(2*wfsn.nValidLenslet) - tt*pinv(tt)); %NGS slopes -> NGS tip-tilt (72x72)
            

%10\ PHASE RECONSTRUCTOR
bif    = gaussianInfluenceFunction(0.2,d); % Gaussian influence function
dm   = deformableMirror(nL+1,'modes',bif,'resolution',tel.resolution); % dm class
calib = calibration(dm,wfs,sref,sref.wavelength/40,1); % Interaction matrix/command matrix calibration
calib.cond = 30; % conditionning of the MI inversion
mc = calib.M; % the reconstructor slopes->command
nActu = size(mc,1);
Fmat = (tel.pupil(:).*dm.modes.modes)*mc; % the phase reconstructor slopes->phase


%% RUN THE END-TO-END SIMULATION

% simulation parameters
dt    = 0;
count = 0;
nExp = 1;%simulation time in seconds.
nIter = round(nExp*Fe); %simulation time in iterations;

tel = tel + atm; %Associate telescope and atmosphere;
lgs = lgs.*tel*wfs; %Define the light optical path for LGSs;
sref = sref.*tel*ts; %Define the light optical path for the Truth sensor;
ngs = ngs.*tel*wfsn; %Define the light optical path for the Truth sensor;

% Instantiate outputs
% Slopes time-vector
tsSl = zeros(nSl,nIter);
wfsSl = zeros(nSl,nLgs,nIter);
wfsnSl = zeros(nSl,nNgs,nIter);
tomoSl_LA = zeros(nSl,nIter);
tomoSl_LGS = zeros(nSl,nIter);
tomoSl_GLAO = zeros(nSl,nIter);
% WFSs camera
tsCam = zeros(nPx*nL,nPx*nL,nIter);
wfsCam = zeros(nPx*nL,nL*nPx*nLgs,nIter);
wfsnCam = zeros(nPx_n*nL_n,nPx_n*nL_n*nNgs,nIter);
% Point-wise phase
phLgs  = zeros(resTel,resTel,nLgs,nIter);
phTs  = zeros(resTel,resTel,nIter);
phNgs  = zeros(resTel,resTel,nNgs,nIter);
% OPD
wfeLgs = zeros(nLgs,nIter);
wfeTs   = zeros(1,nIter);
wfeNgs = zeros(nNgs,nIter);
wfeTomo_LA   = zeros(1,nIter);
wfeRes_LA   = zeros(1,nIter);
wfeTomo_LGS   = zeros(1,nIter);
wfeRes_LGS   = zeros(1,nIter);
wfeTomo_GLAO   = zeros(1,nIter);
wfeRes_GLAO   = zeros(1,nIter);

display = false;
% if true 1.5s/iteration
% if false 0.4s/iteration

close all;
if exist('hwait')
    hwait.delete;
end
hwait = waitbar(0,'Loop is being closed...');

for kIter=1:nIter
    tic;
    %1\ Updating phase screens
    +tel;
    
     %2\ Propagating sources to the TS
    +sref;
    % store phase
    phTs(:,:,kIter) = sref.meanRmPhase;
    wfeTs(kIter)    =  std(ts.slopes*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    %store TS slopes
    tsSl(:,kIter) = ts.slopes;
    tsCam(:,:,kIter) = ts.camera.frame;
    
            
    %3\ Propagating sources to the LGS WFS
    +lgs;    
    % store phase at LGSs positions
    for j=1:nLgs
        phLgs(:,:,j,kIter) = lgs(j).meanRmPhase;
        wfeLgs(j,kIter)    =  std(wfs.slopes(:,j)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    end    
    %store slopes
    wfsSl(:,:,kIter) = wfs.slopes;
    wfsCam(:,:,kIter) = wfs.camera.frame;
    
    
    %4\ Propagating sources to the NGS WFS
    +ngs;    
    % store phase at LGSs positions
    for j=1:nNgs
        phNgs(:,:,j,kIter) = ngs(j).meanRmPhase;
        wfeNgs(j,kIter)    =  std(wfsn.slopes(:,j)*wfsnPscale_simu*constants.arcsec2radian*D*1e9/wfsn.slopesUnits);
    end    
    %store slopes
    wfsnSl(:,:,kIter) = wfsn.slopes;
    wfsnCam(:,:,kIter) = wfsn.camera.frame;
    
    
    %5\ Joint MMSE tomography
    tomoSl_LA(:,kIter) = R_LA*[wfsn.slopes(:);wfs.slopes(:)];   
    wfeTomo_LA(kIter)    =  std(tomoSl_LA(:,kIter)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    wfeRes_LA(kIter)    =  std((tomoSl_LA(:,kIter)-ts.slopes)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
           
    %6\ Split MMSE tomography
    tomoSl_LGS(:,kIter) = R_LGS*[wfs.slopes(:)] + R_TT*mean(wfsn.slopes,2);   
    wfeTomo_LGS(kIter)    =  std(tomoSl_LGS(:,kIter)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    wfeRes_LGS(kIter)    =  std((tomoSl_LGS(:,kIter)-ts.slopes)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    
    %6\ GLAO
    tomoSl_GLAO(:,kIter) = mean(wfs.slopes,2) + R_TT*mean(wfsn.slopes,2);   
    wfeTomo_GLAO(kIter)    =  std(tomoSl_GLAO(:,kIter)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    wfeRes_GLAO(kIter)    =  std((tomoSl_GLAO(:,kIter)-ts.slopes)*wfsPscale_simu*constants.arcsec2radian*d*1e9/wfs.slopesUnits);
    
    if display
        displayResults(42,phLgs,phTs,wfsCam,tsCam,wfeLgs,wfeTs,wfeTomo_LA,wfeRes_LA,kIter,nIter,nExp,nPxWfs)
    end
    
    
    % Updating the waiting bar and remaining time
    waitbar(kIter/nIter);
    dt  = dt + toc();
    if kIter>1 && kIter<nIter
        fprintf(1, repmat('\b',1,count)); %delete line before
        count = fprintf('\nRemaining simulation time: %0.5g s\n',dt*(nIter-kIter)/kIter);
    elseif kIter==nIter
        fprintf(1, repmat('\b',1,count)); %delete line before
        count = fprintf('\nRemaining simulation time: %0.5g s\n',dt*(nIter-kIter)/kIter);
    end
end
close(hwait);

if ~display
    displayResults(42,phLgs,phTs,wfsCam,tsCam,wfeLgs,wfeTs,wfeTomo_LA,wfeRes_LA,kIter,nIter,nExp,nPxWfs)
    %display results at the end of the simulations
end
%% Compare reconstructed phase

%1\ Reconstructing phase maps
phTs = reshape(Fmat*tsSl,resTel,resTel,nIter); % TS phase
phTomo_LA = reshape(Fmat*tomoSl_LA,resTel,resTel,nIter); % Joint Tomography phase
phTomo_LGS = reshape(Fmat*tomoSl_LGS,resTel,resTel,nIter); % Split Tomography phase
phTomo_GLAO = reshape(Fmat*tomoSl_GLAO,resTel,resTel,nIter); % GLAO phase

vTs = zeros(1,nIter);
vTomo_LA = zeros(1,nIter);
vTomo_LGS = zeros(1,nIter);
vTomo_GLAO = zeros(1,nIter);
vRes_LA = zeros(1,nIter);
vRes_LGS = zeros(1,nIter);
vRes_GLAO = zeros(1,nIter);

close all;
figure
for k = 1:nIter
    subplot(3,4,1)
    imagesc(squeeze(phTs(:,:,k)))
    title('Truth sensor','interpreter','latex','FontSize',18);
    set(gca,'XTick',[],'Ytick',[]);
    pbaspect([1,1,1]);
    subplot(3,4,2)
    imagesc(squeeze(phTomo_LA(:,:,k)))
    title('Joint tomography','interpreter','latex','FontSize',18);
    set(gca,'XTick',[],'Ytick',[]);
    pbaspect([1,1,1]);
     subplot(3,4,3)
    imagesc(squeeze(phTomo_LGS(:,:,k)))
    title('Split tomography','interpreter','latex','FontSize',18);
    pbaspect([1,1,1]);
    set(gca,'XTick',[],'Ytick',[]);
     subplot(3,4,4)
    imagesc(squeeze(phTomo_GLAO(:,:,k)))
    title('GLAO','interpreter','latex','FontSize',18);
    pbaspect([1,1,1]);
    set(gca,'XTick',[],'Ytick',[]);
    
     subplot(3,4,6)
    imagesc(squeeze(phTs(:,:,k) - phTomo_LA(:,:,k)))
    title('Residual','interpreter','latex','FontSize',18);
    pbaspect([1,1,1]);
    set(gca,'XTick',[],'Ytick',[]);
     subplot(3,4,7)
    imagesc(squeeze(phTs(:,:,k) - phTomo_LGS(:,:,k)))
    title('Residual','interpreter','latex','FontSize',18);
    pbaspect([1,1,1]);
    set(gca,'XTick',[],'Ytick',[]);
    subplot(3,4,8)
    imagesc(squeeze(phTs(:,:,k) - phTomo_GLAO(:,:,k)))
    title('Residual','interpreter','latex','FontSize',18);
    pbaspect([1,1,1]);
    set(gca,'XTick',[],'Ytick',[]);
    
    vTs(k) = std(reshape(phTs(:,:,k),resTel^2,1))*2e9;
    vTomo_LA(k) = std(reshape(phTomo_LA(:,:,k),resTel^2,1))*2e9;
    vRes_LA(k) = std(reshape(phTs(:,:,k)-phTomo_LA(:,:,k),resTel^2,1))*2e9;
    vTomo_LGS(k) = std(reshape(phTomo_LGS(:,:,k),resTel^2,1))*2e9;
    vRes_LGS(k) = std(reshape(phTs(:,:,k)-phTomo_LGS(:,:,k),resTel^2,1))*2e9;
    vTomo_GLAO(k) = std(reshape(phTomo_GLAO(:,:,k),resTel^2,1))*2e9;
    vRes_GLAO(k) = std(reshape(phTs(:,:,k)-phTomo_GLAO(:,:,k),resTel^2,1))*2e9;
    
    subplot(3,3,7:9)
    plot((1:k)*nExp/nIter,vTs(1:k),'k');hold on;
    plot((1:k)*nExp/nIter,vTomo_LA(1:k),'b');
    plot((1:k)*nExp/nIter,vTomo_LGS(1:k),'r');
    plot((1:k)*nExp/nIter,vTomo_GLAO(1:k),'m');
    plot((1:k)*nExp/nIter,vRes_LA(1:k),'b--');
    plot((1:k)*nExp/nIter,vRes_LGS(1:k),'r--');
    plot((1:k)*nExp/nIter,vRes_GLAO(1:k),'m--');
    xlabel('Ellapsed time (s)','interpreter','latex','FontSize',18);
    ylabel('OPD (nm)','interpreter','latex','FontSize',18);
    legend({'Truth sensor','Joint Tomography','Split Tomography','GLAO','Joint residual','Split residual','GLAO residual'},'interpreter','latex','FontSize',10,'Location','southwest');
    ylim([0,500]);
    drawnow;
    pause(0.05)
end

