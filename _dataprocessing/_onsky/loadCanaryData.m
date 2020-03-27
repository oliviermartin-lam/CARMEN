%% MANAGE WORKING PATHS
clear all;
close all;
clc;
path_root = '/run/media/omartin/DATADISK/DATA/CANARY_DATA/'; %to be changed
path_work = '/home/omartin/Projects/CARMEN/CODES'; %to be changed
addpath(path_root,genpath(path_work));


%% MANAGE DATA PATHS
date         = '20130917'; % to be changed
path_data    = [path_root,date];
path_im      = [path_root,date];
path_calib   = [path_root,'CALIBRATION/'];
path_ir      = [path_data,'/ir/'];     

%% GRAB DATA ID
d = dir(fullfile(path_ir, '*script274*'));
im_id = fullfile(path_ir, {d.name});
nObj   = length(im_id);

%% RESTORE AND PROCESS THE IMAGE
% SCAO: 1-4-7-10 -13
% GLAO: 2-5-8-11 -14
% MOAO: 3-6-9-12 -15
k_id = 2; %object number
tmp = im_id{k_id};
ll = length(path_ir);
obj_name = tmp(ll+15:ll+23);       
trs = handleCanaryData(date,path_data,path_im,path_calib,'obj_name',obj_name);
