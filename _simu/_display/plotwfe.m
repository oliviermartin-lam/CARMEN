clear all;close all;clc;
path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles_old/';
path_carmen =  '/home/omartin/Projects/CARMEN/';
addpath(genpath(path_oomao),genpath(path_carmen));
path_res = '/home/omartin/Projects/CARMEN/RES/CNN_CANARY/';

%% LOAD fits filetomo_MMSE
snoise = '';
test = 3;

switch test
    case 1
        path_test = 'Test_1/';
        r0 = 16;
    case 2
        path_test = 'Test_2/';
        r0 = 12;
    case 3
        path_test = 'Test_3/';
        r0 = 8;
end


slopes_MMSE = fitsread([path_res,path_test,'Test_r0_',num2str(r0),'_MMSE_',snoise,'noise.fits']);
slopes_MLP = fitsread([path_res,path_test,'Test_r0_',num2str(r0),'_MLP_',snoise,'noise.fits'])';
slopes_Conv = fitsread([path_res,path_test,'Test_r0_',num2str(r0),'_Conv_',snoise,'noise.fits'])';
slopes_TS = fitsread([path_res,path_test,'Test_r0_',num2str(r0),'_Ideal_',snoise,'noise.fits'])';

slopes_MMSE = bsxfun(@minus,slopes_MMSE,mean(slopes_MMSE,2));
slopes_MLP =  bsxfun(@minus,slopes_MLP,mean(slopes_MLP,2));
slopes_Conv =  bsxfun(@minus,slopes_Conv,mean(slopes_Conv,2));
slopes_TS =  bsxfun(@minus,slopes_TS,mean(slopes_TS,2));

nExp = size(slopes_TS,2);

% Get the tomographic erors error
tomo_MMSE = slopes_TS - slopes_MMSE;
tomo_MLP = slopes_TS - slopes_MLP;
tomo_CNN = slopes_TS - slopes_Conv;

% get wavefront error
wfe = [getWaveFrontErrorFromSlopes_CANARY(tomo_MMSE)...
    getWaveFrontErrorFromSlopes_CANARY(tomo_MLP)...
    getWaveFrontErrorFromSlopes_CANARY(tomo_CNN)];

% Get the covariance matrix
C_MMSE = tomo_MMSE*tomo_MMSE'/nExp;
C_MLP = tomo_MLP*tomo_MLP'/nExp;
C_CONV = tomo_CNN*tomo_CNN'/nExp;

clear slopes_TS slopes_MMSE slopes_Conv slopesMLP tomo_MMSE tomo_MLP tomo_Conv

% Create a phase reconstructor
 parFileCanary_3NGS
% % telescope, source
 tel = telescope(D,'resolution',resTel,'obstructionRatio',cobs);

wvl = 1.65e-6;
waveNum = 2*pi/wvl;
psTS = 0.219245;
R =  waveNum*(0.5 * tel.D  * psTS / (3600.*180./pi)) * fitsread('GLOB_mrz_7x7.fits')'; %in nm
C_MMSE = R*C_MMSE*R';
C_MLP = R*C_MLP*R';
C_CONV = R*C_CONV*R';
z = zernike(2:size(R ,1)+1,'resolution',tel.resolution);

% Get the PSF
fov = 512;
Samp = 206264.8*1e3*wvl/tel.D/30/2;
psInMas =constants.radian2mas*wvl/tel.D/Samp/2;

psf_MMSE = tools.otf2psf(tools.modes2Otf(C_MMSE,z.modes,tel.pupil,fov,Samp));
psf_MLP = tools.otf2psf(tools.modes2Otf(C_MLP,z.modes,tel.pupil,fov,Samp));
psf_CONV = tools.otf2psf(tools.modes2Otf(C_CONV,z.modes,tel.pupil,fov,Samp));
psf_DL = tools.otf2psf(tools.modes2Otf(0*C_CONV,z.modes,tel.pupil,fov,Samp));
S = max(psf_DL(:));

%


SR = [tools.getStrehl(psf_MMSE,tel.pupil,Samp) ...
    tools.getStrehl(psf_MLP,tel.pupil,Samp) ...
    tools.getStrehl(psf_CONV,tel.pupil,Samp)];

FWHM = [tools.getFWHM(psf_MMSE,psInMas,4,'contour') ...
    tools.getFWHM(psf_MLP,psInMas,4,'contour') ...
    tools.getFWHM(psf_CONV,psInMas,4,'contour')];

% a = tools.getEncircledEnergy(tools.interpolate(psf_DL,2048));
% px = find(a);idx = find(a>0.5);idx = idx(1);
% E50_DL = (px(idx-1) + (0.5 - a(idx-1))/(a(idx) - a(idx-1)))*psInMas;

a = tools.getEncircledEnergy(tools.interpolate(psf_MMSE,2048));
px = find(a);idx = find(a>0.5);idx = idx(1);
E50_MMSE = (px(idx-1) + (0.5 - a(idx-1))/(a(idx) - a(idx-1)))*psInMas;

a = tools.getEncircledEnergy(tools.interpolate(psf_MLP,2048));
px = find(a);idx = find(a>0.5);idx = idx(1);
E50_MLP = (px(idx-1) + (0.5 - a(idx-1))/(a(idx) - a(idx-1)))*psInMas;

a = tools.getEncircledEnergy(tools.interpolate(psf_CONV,2048));
px = find(a);idx = find(a>0.5);idx = idx(1);
E50_CNN = (px(idx-1) + (0.5 - a(idx-1))/(a(idx) - a(idx-1)))*psInMas;

E50 = [E50_MMSE E50_MLP E50_CNN];

[SR;FWHM;E50;wfe]

% DISPLAY
fov = 64;
close all;
figure;
A = abs([tools.crop(psf_MMSE,fov),tools.crop(psf_MLP,fov),tools.crop(psf_CONV,fov)]/S);
s = surf(A);
view(0,0);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
pbaspect([3,1,1]);
cb = colorbar('east');
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 20;
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
set(cb,'position',[0.9 .259 .01 .517])
hold on
SR = max(psf_MMSE(:))/S;
plot3([xlim()],[33,33],[SR,SR],'k--');
SR = max(psf_MLP(:))/S;
plot3([xlim()],[33,33],[SR,SR],'k--');
SR = max(psf_CONV(:))/S;
plot3([xlim()],[33,33],[SR,SR],'k--');


%% HEIGHT TEST
clear all;close all;clc;
path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles_old/';
path_carmen =  '/home/omartin/Projects/CARMEN/';
addpath(genpath(path_oomao),genpath(path_carmen));
path_res = '/home/omartin/Projects/CARMEN/RES/CNN_CANARY/';

wfe_cnn = zeros(1,15);
wfe_mlp = zeros(1,15);
wfe_mmse = zeros(1,15);
for h = 1:15
    % Get fits file
    CNN_h = fitsread([path_res,'Test_heights/Conv_height_',num2str(h-1),'.fits'])';
    CNN_h = CNN_h - mean(CNN_h,2);
    MLP_h = fitsread([path_res,'Test_heights/MLP_height_',num2str(h-1),'.fits'])';
    MLP_h = MLP_h - mean(MLP_h,2);
    MMSE_h = fitsread([path_res,'Test_heights/tomo_slopes_',num2str(h*1000),'km.fits']);
    MMSE_h = MMSE_h - mean(MMSE_h,2);
    Truth_h = fitsread([path_res,'Test_heights/Ideal_height_',num2str(h-1),'.fits'])';
    Truth_h = Truth_h - mean(Truth_h,2);
    % Get the wfe
    wfe_mmse(h) = getWaveFrontErrorFromSlopes_CANARY(MMSE_h - Truth_h);
    wfe_mlp(h) = getWaveFrontErrorFromSlopes_CANARY(MLP_h - Truth_h);
    wfe_cnn(h) = getWaveFrontErrorFromSlopes_CANARY(CNN_h - Truth_h);
end

wfe_mmse2 = fitsread([path_res,'Test_heights/mmse_perf.fits']);

%% WFE = f(h)
h = 1:1:15;
pp = polyfit(h,wfe_mmse,1);
pp2 = polyfit(h,wfe_mmse2(4,:),1);
wfe_th = polyval(pp,h);

close all;
figure;
plot(h,wfe_mlp,'bo--','MarkerFaceColor','b','MarkerSize',7);hold on
plot(h,wfe_cnn,'rd--','MarkerFaceColor','r','MarkerSize',7);
%plot(h,wfe_mmse2(1,:) + ( wfe_th(2) - wfe_mmse2(1,2)),'ks--','MarkerFaceColor','k','MarkerSize',7);
plot(h,wfe_mmse2(2,:) + ( wfe_th(7) - wfe_mmse2(2,7)),'md--','MarkerFaceColor','m','MarkerSize',7);
plot(h,wfe_mmse2(3,:) + ( wfe_th(12) - wfe_mmse2(3,12)),'ko--','MarkerFaceColor','k','MarkerSize',7);

plot(h,wfe_th,'k--','linewidth',1.5);
xlabel('Altitude layer (km)','interpreter','latex','fontsize',24);
ylabel('Wave front error (nm)','interpreter','latex','fontsize',24);
legend({'MLP','CNN','MMSE 7km','MMSE 12km','Theoretical limit'},...
    'interpreter','latex','fontsize',20,'Location','northwest');
pbaspect([1,1,1]);
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');


%% WFE = f(h) - CNN vs MMSE
h = 1:1:15;
pp = polyfit(h,wfe_mmse,1);
pp2 = polyfit(h,wfe_mmse2(4,:),1);
wfe_th = polyval(pp,h);

close all;
figure;
plot(h,wfe_cnn,'bo--','MarkerFaceColor','b','MarkerSize',7);hold on
plot(h,wfe_mmse2(2,:) + ( wfe_th(7) - wfe_mmse2(2,7)),'rs--','MarkerFaceColor','r','MarkerSize',7);
plot(h,wfe_th,'k--','linewidth',1.5);
xlabel('Altitude layer (km)','interpreter','latex','fontsize',24);
ylabel('Wave front error (nm)','interpreter','latex','fontsize',24);
legend({'CNN','MMSE','Theoretical limit'},...
    'interpreter','latex','fontsize',20,'Location','northwest');
pbaspect([1,1,1]);
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
