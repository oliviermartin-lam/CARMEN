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
path_simu       = [path_carmen,'_simu/']; % the simulation folder path
path_processing = [path_carmen,'_dataprocessing/']; % the simulation folder path
path_scidar     = [path_carmen,'_scidar/']; % the directory containing Scidar data
addpath(genpath(path_oomao),genpath(path_simu),genpath(path_processing),path_scidar);

%% MANAGE OPTIONS
dataType  = 'test'; 
% Can be : 

% training : 2-layers atmosphere with a variable-altitude for the second layer. The altitude is linearly dsampled from hmin to hmax using nScreens values (parFile).
% Do nIter independent simulations (defined in parFile) of slopes

% test : user-defined atmosphere (defined in parFile). Simulate time-series  (frozen-flow turbulence) of slopes for a specific atmosphere

% obs : Stereo-scidar-based simulations. The user must define the temporal
% resolution and the number of reconstructed layers (parFile). Simulate nBins time-series of slopes

flagNoise = false; %% TBC %%
flagSave  = false;
simuCase  = 'Canary_4LGS'; % TBC
flagDisp  = false; % if true, some figures will pop up
frozenflow= false; % if true, the code simulates temporally correlated time-series of slopes accounting for the frozen-flow assumption
getZernike= false; % if true, reconstruct Zernike coefficients from phase/slopes
flagMMSE  = false; %% TBC %%

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
    gsRef   = [];
else
    % LASER GUIDE STAR
    nPhoton = mean(photoGs.zeroPoint*10.^(-0.4*magGs));
    if nLayerNa > 0
        hNa     = linspace(hGs-naThickness/2,hGs+naThickness/2,nLayerNa);
    else
        naWeight= 1;
        hNa     = hGs;
    end
    gs      = laserGuideStar(tel.D/nL,tel.D,hGs,spotFwhm,nPhoton,naWeight/sum(naWeight),'zenith',rGs,'azimuth',dGs,'height',...
        hNa,'wavelength',photoGs,'viewPoint',viewPoint);
    gsRef   = source('magnitude',magGs,'zenith',rGs,'azimuth',dGs,'wavelength',photoGs);
    
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
if flagDisp
    figure;
    pp = gs.polar;
    for i=1:nGs
        pp(i).MarkerFaceColor = 'b';
        pp(i).MarkerSize= 10;
        pp(i).Marker= 'o';
        pp(i).Color   = 'b';
    end
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

% FROM SHACK-HARTMAN WFS SLOPES/PHASE TO ZERNIKE COEFFICIENTS
[S2Z,P2Z]   = calibrateZernikeReconstructionMatrix(tel,ts,sref,'nModes',nZern,'amp',sref.wavelength/40);
close all;

if flagDisp
    % -------------------- plot linearity curves
    jIndex      = [2,4,10,nZern];
    z           = zernike(jIndex,'resolution',tel.resolution,'D',tel.D);
    nTest       = 20;
    zMax        = 1e3;
    nZernTest   = numel(jIndex);
    zTrue       = linspace(-zMax,zMax,nTest);
    z_rec1      = zeros(nZernTest,nTest);
    z_rec2      = zeros(nZernTest,nTest);
    for k = 1:nTest
        for j=1:nZernTest
            ph          = reshape(z.modes(:,j)*2*pi*zTrue(k)*1e-9/sref.wavelength,tel.resolution,[]);
            sref        = sref.*tel;
            sref.phase  = ph;
            sref        = sref*ts;
            z_rec1(j,k) = S2Z(jIndex(j)-1,:)* ts.slopes;
            z_rec2(j,k) = P2Z(jIndex(j)-1,:)*ph(:);
        end
    end
    

    figure;
    plot(zTrue,z_rec1(1,:),'ks--','MarkerFaceColor','k','MarkerSize',5);hold on;
    plot(zTrue,z_rec1(2,:),'bo--','MarkerFaceColor','b','MarkerSize',5);
    plot(zTrue,z_rec1(3,:),'rd--','MarkerFaceColor','r','MarkerSize',5);
    plot(zTrue,z_rec1(4,:),'mh--','MarkerFaceColor','m','MarkerSize',5);
    xlabel('Ground truth (nm)','interpreter','latex','fontsize',20);
    ylabel('Slopes-based Reconstruction (nm)','interpreter','latex','fontsize',20);
    set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
    legend({'j=2','j=4','j=10',sprintf('j=%d',nZern)},'interpreter','latex','FontSize',18,'Location','northwest');
    
    figure;
    plot(zTrue,z_rec2(1,:),'ks--','MarkerFaceColor','k','MarkerSize',5);hold on;
    plot(zTrue,z_rec2(2,:),'bo--','MarkerFaceColor','b','MarkerSize',5);
    plot(zTrue,z_rec2(3,:),'rd--','MarkerFaceColor','r','MarkerSize',5);
    plot(zTrue,z_rec2(4,:),'mh--','MarkerFaceColor','m','MarkerSize',5);
    xlabel('Ground truth (nm)','interpreter','latex','fontsize',20);
    ylabel('Phase-based Reconstruction (nm)','interpreter','latex','fontsize',20);
    set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
    legend({'j=2','j=4','j=10',sprintf('j=%d',nZern)},'interpreter','latex','FontSize',18,'Location','northwest');
    
end

%% CREATING SAVING FOLDERS
if flagSave
    st = {'WFS_CAM','WFS_SLOPES','TS_CAM','TS_SLOPES','TOMO_SLOPES','ZERNIKE_COEFFICIENTS'};
    saveDir = [path_res,upper(simuCase),'/',dataType,'/'];
    for k=1:numel(st)
        path = [saveDir,st{k}];
        if ~isfolder(path)
            mkdir(path)
        end
    end
end
%% RUN THE END-TO-END SIMULATION


switch dataType
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'training'
        alt = linspace(hmin/1e3,hmax/1e3,nScreens);
        for l=1:nScreens
            atm = atmosphere(photoGs,r0,mean(L0_training),'layeredL0',L0_training,'fractionnalR0',fractionalR0_training,...
                'altitude',[0,alt(l)*1e3],'windSpeed',10*ones(1,2),'windDirection',zeros(1,2));
            t0  = tic();
            trs = generateTelemetry(tel,atm,gs,sref,wfs,ts,nIter,'frozenflow',false,'ron',ron,...
                'mmse',flagMMSE,'S2Z',S2Z,'P2Z',P2Z,'getZernike',getZernike,'gsRef',gsRef);
            tf  = toc(t0);
            fprintf('Done in %.1f s',tf)
            
            
            if flagSave
                fitswrite(trs.wfsSl,[saveDir,'WFS_SLOPES/offaxiswfss_slopes_',num2str(atm.layer(2).altitude),'km_noise_',num2str(flagNoise),'.fits']);
                fitswrite(trs.tsSl,[saveDir,'TS_SLOPES/ts_slopes_',num2str(atm.layer(2).altitude),'km_noise_',num2str(flagNoise),'.fits']);
                fitswrite(trs.wfsCam,[saveDir,'WFS_CAM/offaxiswfss_cam_',num2str(atm.layer(2).altitude),'km_noise_',num2str(flagNoise),'.fits']);
                fitswrite(trs.tsCam,[saveDir,'TS_CAM/ts_slopes_',num2str(atm.layer(2).altitude),'km_noise_',num2str(flagNoise),'.fits']);
                if flagMMSE
                    fitswrite(trs.tomoSl,[saveDir,'TOMO_SLOPES/tomo_slopes_',num2str(atm.layer(2).altitude),'km_noise_',num2str(flagNoise),'.fits']);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'test'
        
        % GENERATING TEST TIME_SERIES FOR A SPECIFIC ATMOSPHERE
        t0  = tic();
        trs = generateTelemetry(tel,atm,gs,sref,wfs,ts,10,'frozenflow',false,'ron',ron,...
            'mmse',flagMMSE,'S2Z',S2Z,'P2Z',P2Z,'getZernike',getZernike,'gsRef',gsRef);
        tf  = toc(t0);
        fprintf('Done in %.1f s\n',tf)
        
        if ~exist('atmName','var')
            atmName = 'noname';
        end
        if flagSave
            fitswrite(trs.wfsSl,[saveDir,'WFS_SLOPES/offaxiswfss_slopes_',atmName,'_noise_',num2str(flagNoise),'.fits']);
            fitswrite(trs.tsSl,[saveDir,'TS_SLOPES/ts_slopes_',atmName,'_noise_',num2str(flagNoise),'.fits']);
            fitswrite(trs.wfsCam,[saveDir,'WFS_CAM/offaxiswfss_cam_',atmName,'_noise_',num2str(flagNoise),'.fits']);
            fitswrite(trs.tsCam,[saveDir,'TS_CAM/ts_slopes_',atmName,'_noise_',num2str(flagNoise),'.fits']);
            if flagMMSE
                fitswrite(trs.tomoSl,[saveDir,'TOMO_SLOPES/tomo_slopes_',num2str(atm.nLayer),'layers_noise_',num2str(flagNoise),'.fits']);
            end
        end
        
        if isfield(trs,'wfe')
            trs.wfe.uncorrected
            trs.wfe.uncorrected_th
            if flagMMSE
                trs.wfe.mmse
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'obs'
        % GET CN2 PROFILES
        [Cn2_c,alt_c,r0,date_bin]  = readScidarData(path_scidar,profileConfig,dt,nL_c,photoGs.wavelength,flagDisp);
        nBin   = size(Cn2_c,2);
        vmean  = mean(windSpeed);
        
        % SIMULATING TIME-SERIES FOR EACH PROFILE
        for kBin = 1:nBin
            % UPDATING ATMOSPHERE
            idGood = Cn2_c(kBin,:) ~=0;
            wK     = Cn2_c(kBin,idGood);
            wK     = wK/sum(wK);
            atm   = atmosphere(photoGs,r0(kBin),mean(L0),'layeredL0',mean(L0),'fractionnalR0',wK,...
                'altitude',alt_c(kBin,idGood),'windSpeed',vmean*ones(1,nL_c),'windDirection',zeros(1,nL_c));
            
            % GENERATING TELEMETRY
            t0  = tic();
            trs = generateTelemetry(tel,atm,gs,sref,wfs,ts,nIter,'frozenflow',false,'ron',ron,...
                'mmse',false,'S2Z',S2Z,'P2Z',P2Z,'getZernike',getZernike,'gsRef',gsRef);
            tf  = toc(t0);
            fprintf('Done in %.1f s\n',tf)
            
            %MMSE RECONSTRUCTOR
            if flagMMSE
                if kBin == 1
                    Rmmse = getMMSE(tel,atm,gs,wfs,sref,S2Z);
                end
                trs.tomoSl = Rmmse*reshape(trs.wfsSl,nSl*nGs,nIter);
            end
            
            % SAVE TELEMETRY
            if flagSave
                fitswrite(trs.wfsSl,[saveDir,'WFS_SLOPES/''offaxiswfss_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                fitswrite(trs.tsSl,[saveDir,'TS_SLOPES/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                fitswrite(trs.wfsCam,[saveDir,'WFS_CAM/offaxiswfss_cam_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                fitswrite(trs.tsCam,[saveDir,'TS_CAM/ts_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                if flagMMSE
                    fitswrite(trs.tomoSl,[saveDir,'TOMO_SLOPES/tomo_slopes_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                end
                if getZernike
                    zerAll = zeros(2*nGs+2,nZern,nIter);
                    zerAll(1,:,:)     = trs.zer.ts_truth;
                    zerAll(2,:,:)     = trs.zer.ts_rec;
                    zerAll(3:2+nGs,:,:) = trs.zer.wfs_truth;
                    zerAll(3+nGs:end,:,:) = trs.zer.wfs_rec;
                    fitswrite(zerAll,[saveDir,'ZERNIKE_COEFFICIENTS/zer_coeffs_',num2str(nL_c),'bins_',num2str(dt),'mn_',profileConfig,'_',num2str(date_bin(kBin)),'.fits']);
                end
            end
        end
        
    otherwise
        error('The data type to be simulated is not understood')
end


%%
if flagDisp
    close all;clc;
    fontsize = 20;
    
    figure;
    imagesc(trs.wfsCam(:,1:1:nL*nPx,end)) ;hold on;
    for j=1:nL-1
        plot([xlim()],j*nPx*[1,1],'w');
        plot(j*nPx*[1,1],[ylim()],'w');
    end
    th = 0:pi/50:2*pi;
    x = nL*nPx * cos(th)/2 + nL*nPx/2;
    y = nL*nPx * sin(th)/2 + nL*nPx/2;
    xo = tel.obstructionRatio*nL*nPx * cos(th)/2 + nL*nPx/2;
    yo = tel.obstructionRatio*nL*nPx * sin(th)/2 + nL*nPx/2;
    plot(x,y,'w--');
    plot(xo,yo,'w--');
    pbaspect([1,1,1]);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    xlabel('Detector pixel','interpreter','latex','fontsize',fontsize)
    ylabel('Detector pixel','interpreter','latex','fontsize',fontsize)
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontsize;
    
    figure
    imagesc(1e9*trs.opdNGS(:,:,1,end));hold on;
    pbaspect([1,1,1]);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    xlabel('Detector pixel','interpreter','latex','fontsize',fontsize)
    ylabel('Detector pixel','interpreter','latex','fontsize',fontsize)
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontsize;
    
    
    tmp = trs.wfsSl(1:36,1,end);
    sx = zeros(nL);
    sx(wfs.validLenslet) = tmp;
    tmp = trs.wfsSl(37:72,1,end);
    sy = zeros(nL);
    sy(wfs.validLenslet) = tmp;
    figure
    imagesc([sx,sy]);hold on;
    pbaspect([2,1,1]);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    xlabel('WFS subaperture','interpreter','latex','fontsize',fontsize)
    ylabel('WFS subaperture','interpreter','latex','fontsize',fontsize)
    xticks(1:14)
    xticklabels({'1' '2' '3'  '4'  '5' '6' '7' '1' '2' '3' '4' '5' '6' '7'})
    ax = gca;
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontsize;
end