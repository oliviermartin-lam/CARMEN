clear all;close all;clc;
path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles_old/';
path_carmen =  '/home/omartin/Projects/CARMEN/';
addpath(genpath(path_oomao),genpath(path_carmen));

%%
fov = 32;
path_root = '/run/media/omartin/OlivierMartinHDD/DATA/CANARY_DATA/'; %to be changed
path_calib   = [path_root,'CALIBRATION/'];
badPixMap = fitsread([path_calib,'locateDeadPixelsCAMICAZ.fits']);
samp = constants.radian2mas*1.65e-6/4.2/30/2;
ir = fitsread([path_root,'/20130722/ir/ir_2013-07-23_01h50m27s_ANN.fits']);
irbg =fitsread([path_root,'/20130722/irbg/irbg_2013-07-23_00h11m12s_test.fits']);
psf_ANN = tools.processImage(ir,irbg,1,badPixMap,samp,'fov',fov,'rebin',4);
psf_ANN = psf_ANN/sum(psf_ANN(:));

ir = fitsread([path_root,'/20130722/ir/ir_2013-07-23_01h48m01s_MOAO.fits']);
psf_MOAO = tools.processImage(ir,irbg,1,badPixMap,samp,'fov',fov,'rebin',4);
psf_MOAO = psf_MOAO/sum(psf_MOAO(:));

ir = fitsread([path_root,'/20130722/ir/ir_2013-07-22_17h41m59s_FP_psf_reference.fits']);
irbg =fitsread([path_root,'/20130722/irbg/irbg_2013-07-22_17h40m52s_FP_psf_reference.fits']);
psf_NCPA = tools.processImage(ir,irbg,1,badPixMap,samp,'fov',fov,'rebin',4);

psf_NCPA = psf_NCPA / sum(psf_NCPA(:));
S = max(psf_NCPA(:));
close all;
figure;
A = 1e2*abs([psf_MOAO,psf_ANN]/S);
s = surf(A);
view(0,0);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
pbaspect([3,1,1]);
cb = colorbar('east');
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 20;
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
set(cb,'position',[0.9 .259 .01 .517])


[tools.getFWHM(psf_MOAO,psInMas,4,'contour') ...
    tools.getFWHM(psf_ANN,psInMas,4,'contour') ...
  ]
