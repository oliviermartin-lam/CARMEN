%% MANAGE WORKING PATHS
clear all;
close all;
clc;
path_root = '/run/media/omartin/OlivierMartinHDD/DATA/CANARY_DATA/'; %to be changed
path_work = '/home/omartin/Projects/CARMEN/'; %to be changed
addpath(genpath(path_work));


%% GRAB DATA FROM SCRIPT ID

% MANAGE DATA PATHS
date         = '20130917'; % to be changed
path_data    = [path_root,date];
path_calib   = [path_root,'CALIBRATION/'];
path_ir      = [path_data,'/ir/'];     

% Get fits file names
d = dir(fullfile(path_ir, '*script274*'));
im_id = fullfile(path_ir, {d.name});
nObj   = length(im_id);

% RESTORE AND PROCESS THE IMAGE
% SCAO: 1-4-7-10 -13
% GLAO: 2-5-8-11 -14
% MOAO: 3-6-9-12 -15

% select individual file
k_id = 2; %object number
tmp = im_id{k_id};
ll = length(path_ir);
obj_name = tmp(ll+15:ll+23);       
trs = handleCanaryData(date,path_data,path_calib,'obj_name',obj_name);

%% GRAB DATA FROM FITS NAME

date         = '20130723'; % to be changed
path_data    = [path_root,date];
path_calib   = [path_root,'CALIBRATION/'];
path_ir      = [path_data,'/ir/'];     
nBox         = 64;

trs = handleCanaryData(date,path_data,path_calib,'obj_name','22h34m59s','getImageOnly',true,'resolution',nBox);
psf_ol = trs.im_sky;
trs = handleCanaryData(date,path_data,path_calib,'obj_name','22h35m25s','getImageOnly',true,'resolution',nBox);
psf_cl = trs.im_sky;

% display
A = [psf_ol,psf_cl];

close all;
figure;
imagesc(log10(A.*(A>0)),[1,4.0])
pbaspect([2,1,1])
set(gca,'XTick',[],'YTick',[]);
cb = colorbar();
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 20;


text(nBox/3.5,nBox/10,'W/O Adaptive Optics','color','white','interpreter','latex','fontsize',20)
text(4.8*nBox/3.5,nBox/10,'With Adaptive Optics','color','white','interpreter','latex','fontsize',20)
text(nBox*2 - 25,nBox-1,'Credit: O. Beltramo-Martin','color','white','interpreter','latex','fontsize',16)


xx = linspace(0,nBox/2,nBox/2)*trs.psInMas;
figure
semilogy(xx,radial(psf_ol),'b');
hold on
semilogy(xx,radial(psf_cl),'r');
pbaspect([1,1,1])
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex' );
xlabel('Angular separation (mas)','interpreter','latex','fontsize',20);
ylabel('Azimuthal profile','interpreter','latex','fontsize',20);
legend({'W/O Adaptive Optics','With Adaptive Optics'},'interpreter','latex','FontSize',20,'Location','northeast');