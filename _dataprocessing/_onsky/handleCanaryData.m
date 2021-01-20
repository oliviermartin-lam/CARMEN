    classdef handleCanaryData < handle
    
    properties  (SetObservable=true)
        % ---------------------------- PATHS ---------------------------- %
        date;
        obj_name;
        path_data;
        path_calib;
        path_profile;
        % -------------------------- TELEMETRY -------------------------- %
        aomode;
        slopes;
        tipTilt;
        hodm_pos;
        ittm_pos;
        focus;
        MI;
        R;
        Rtt;
        Rtomo;
        waveFront;
        % ------------------------ SCIENCE IMAGE ------------------------ %
        im_sky;
        resolution;
        airmass=1;
        wvl;
        psInMas;
        Samp;
        hdr;
        % ------------------------- CALIBRATION ------------------------- %
        staticMap;
        dmIF;
        psf_ncpa;
        otf_ncpa;
        % ------------------------- LOOP STATUS ------------------------- %
        holoop_gain
        holoop_freq ;
        holoop_lat;
        holoop_rtf;
        holoop_ntf;
        holoop_pn;
        ttloop_gain;
        ttloop_freq;
        ttloop_lat;
        ttloop_rtf;
        ttloop_ntf;
        ttloop_pn
    end
    
    
    
    methods
        
        function obj = handleCanaryData(date,path_data,path_calib,varargin)
            inputs = inputParser;
            inputs.addRequired('date', @ischar);
            inputs.addRequired('path_data', @ischar);
            inputs.addRequired('path_calib', @ischar);
            inputs.addParameter('path_profile', [],@ischar);
            inputs.addParameter('resolution',128,@isnumeric);
            inputs.addParameter('obj_name',[],@ischar);
            inputs.addParameter('getImageOnly',false,@islogical);
            inputs.parse(date,path_data,path_calib,varargin{:});
            
            %Inputs
            obj.date         = date;
            obj.path_data    = inputs.Results.path_data;
            obj.path_calib   = inputs.Results.path_calib;
            obj.path_profile = inputs.Results.path_profile;
            obj.resolution   = inputs.Results.resolution;
            obj.obj_name     = inputs.Results.obj_name;
            getImageOnly     = inputs.Results.getImageOnly;
            
            %%
            %GET AND PROCESS THE CAMICAZ IMAGE
            
            %1\ Get the raw image
            path_ir   = [path_data,'/ir/'];
            cd(path_ir);
            path_file= ls(['*',obj.obj_name,'*']);
            path_ir  = [path_ir,path_file];
            path_ir  = path_ir(1:end-1);
            imRaw    = fitsread(path_ir);
            % Get time
            im_header = fitsinfo(path_ir);
            obj.hdr    = im_header.PrimaryData.Keywords;
            list      = obj.hdr(:,1);
            val       = obj.hdr(:,2);
            time_ir  = tools.str2time(obj.obj_name);
            
            %2\ Get the background
            %2.1 Explore files
            path_irbg = [path_data,'/irbg/'];
            cd(path_irbg);
            list     = ls();
            fname    = textscan( list, '%s');
            fname    = fname{1};
            fname =  fname(contains(fname,'.fits'));
            nObj     = length(fname);
            %2.2 Get the time
            time_bg  = zeros(1,nObj);
            for i=1:nObj
                tmp     = fname{i};
                tmp_spt = strsplit(tmp,'_');
                if length(tmp_spt)>1
                    tmp     = tmp_spt{3};
                    time_bg(i) = tools.str2time(tmp);
                end
            end
            [~,closestIndex] = min(abs(time_bg-time_ir));
            %2.3 Load the closest bg file
            ir_bg     = fitsread([path_irbg,fname{closestIndex}]);
            
            %\3 Get the bad pixels map and trs date
            im_header = fitsinfo(path_ir);
            hdr    = im_header.PrimaryData.Keywords;
            list      = hdr(:,1);
            val       = hdr(:,2);
            camName   = upper(cell2mat(val(strcmp(list,'IRCAM'))));
            badPixMap = fitsread([path_calib,'locateDeadPixels',camName,'.fits']);
            
            %4\ Process the image
            obj.wvl     = cell2mat(val(strcmp(list,'FILTER')))*1e-9;
            D           = cell2mat(val(strcmp(list,'TELDIAM')));
            NPIXLD      = cell2mat(val(strcmp(list,'NPIXLD')));
            obj.psInMas = NPIXLD*obj.wvl/D*(1e3*3600*180/pi);
            obj.Samp    = 1/2/NPIXLD;
            obj.airmass   = str2double(cell2mat(val(strcmp(list,'WHTAIRM'))));
            
            badPixMap(imRaw - ir_bg == 0) = 1;
            obj.im_sky = tools.processImage(imRaw,ir_bg,1,badPixMap,obj.Samp,...
                'fovInPixel',obj.resolution,'masking',false,'rebin',0,'tfccd',false,...
                'thresholding',-Inf);
            
            %5\ Observing mode
            obs_mode  = cell2mat(val(strcmp(list,'OBS_MODE')));
            if strcmp(obs_mode,'MOAO')
                rectype = cell2mat(val(strcmp(list,'RECTYPE')));
                if ~isempty(rectype)
                    obj.aomode = upper(rectype);
                else
                    obj.aomode = 'MOAO';
                end
            else
                obj.aomode = "SCAO";
            end
                
                
            if ~getImageOnly
                %5\ Restore the Bench PSF
                %5.1 load the file
                psf_calib = '2013-09-16_21h20m56s_checkPsfQualityBeforeNight';
                time_ncpa = '21h20m56s';
                obj.psf_ncpa  = fitsread([path_data(1:end-8),'20130916/ir/ir_',psf_calib,'.fits']);
                %5.2 load the file
                time_ncpa  = tools.str2time(time_ncpa);
                [~,closestIndex] = min(abs(time_bg-time_ncpa));
                %5.3 Load the closest bg file
                ir_bg = fitsread([path_irbg,fname{closestIndex}]);
                obj.psf_ncpa = tools.processImage(obj.psf_ncpa,ir_bg,1,badPixMap,obj.Samp,...
                    'fovInPixel',obj.resolution,'masking',false,'rebin',0,'tfccd',false,...
                    'thresholding',-Inf);
                %5.4 Get the corresponding OTF with a Nyquist sampling
                tmp = tools.interpolate(obj.psf_ncpa,floor(obj.resolution/obj.Samp));
                tmp = tools.recenterPSF(tmp,4);
                tmp(tmp<0) = 0;
                obj.otf_ncpa = tools.psf2otf(tmp);
                obj.otf_ncpa = obj.otf_ncpa/max(obj.otf_ncpa(:));
                
                %%
                %GET AND PROCESS TELEMETRY
                
                %1\ Get paths
                % look into the ir path
                
                
                %telemetry path
                sl_date   = cell2mat(val(strcmp(list,'SLOPFILE')));
                
                if ~isempty(sl_date)
                    volt_date = cell2mat(val(strcmp(list,'VOLTFILE')));
                    rec_date  = cell2mat(val(strcmp(list,'MC')));
                    mi_date  = cell2mat(val(strcmp(list,'MI')));
                    
                    %2\ Get AO control loop data
                    
                    %2.2 Slopes
                    path_sl   = [path_data,'/slopestl/'];
                    cd(path_sl);
                    path_file  = ls(['*',sl_date,'*']);
                    path_sl    = [path_sl,path_file];
                    path_sl    = path_sl(1:end-1);
                    nSl        = str2num(cell2mat(val(strcmp(list,'WFSNSLO'))));
                    its        = cell2mat(val(strcmp(list,'ITS')));
                    % Prepare reading
                    info   = fitsinfo(path_sl);
                    list_sl= info.Image(1).Keywords(:,1);
                    val_sl = info.Image(1).Keywords(:,2);
                    nWfs   = size(info.Image,2);
                    nExp   = cell2mat(val_sl(strcmp(list_sl,'NAXIS3')));
                    % Prepare compensation of WFS rotation
                    wfs_sym   = str2num(cell2mat(val(strcmp(list,'WFSSYM'))));
                    obj.slopes = zeros(nSl(its),nWfs,nExp);
                    for iWfs=1:nWfs
                        tmp = reshape(fitsread(path_sl,'image',iWfs),nSl(iWfs),[]);
                        obj.slopes(:,iWfs,:)= canaryTools.mirror_SH7(tmp , wfs_sym(iWfs) );
                    end
                    
                    %2.3 DM commandes
                    volt2adu  = 2^15/10.;
                    volt2meter= 0.78*1e-6;
                    adu2volt  = 1./volt2adu;
                    
                    path_volt   = [path_data,'/voltstl/'];
                    cd(path_volt);
                    path_file    = ls(['*',volt_date,'*']);
                    path_volt    = [path_volt,path_file];
                    path_volt    = path_volt(1:end-1);
                    obj.hodm_pos = fitsread(path_volt,'image')';
                    obj.ittm_pos = obj.hodm_pos(end-1:end,:);
                    obj.hodm_pos = obj.hodm_pos(1:end-2,:)*volt2meter;
                    
                    %3\ Get matrices
                    path_rec   = [path_data,'/mc/'];
                    cd(path_rec);
                    path_file  = ls(['*',rec_date,'*']);
                    path_rec   = [path_rec,path_file];
                    path_rec   = path_rec(1:end-1);
                    obj.R      = fitsread(path_rec)'; %in ADU/pixel
                    obj.R      = volt2meter*adu2volt*obj.R;%meter/pixel
                    obj.R      = 0.7804*obj.R(1:end-2,:);
                    % Tip-tilt reconstructor
                    ps_ts      = cell2mat(val(strcmp(list,['PIXARC',num2str(its)])));
                    mrz        = fitsread([path_calib,'/GLOB_mrz_7x7.fits']);
                    obj.Rtt    = mrz(:,1:2)';
                    % Interaction matrix
                    path_mi   = [path_data,'/mi/'];
                    cd(path_mi);
                    path_file  = ls(['*',mi_date,'*']);
                    path_mi   = [path_mi,path_file];
                    path_mi   = path_mi(1:end-1);
                    obj.MI      = fitsread(path_mi)'; %in ADU/pixel
                    obj.MI = obj.MI(:,1:end-2)/volt2meter/0.7804;
                    
                    if ~strcmp(obj.aomode,'scao')
                        % Load the command matrix
                        mct_date    = cell2mat(val(strcmp(list,'MCT')));
                        path_mct    = [path_data,'/mct/'];
                        cd(path_mct);
                        path_file   = ls(['*',mct_date,'*']);
                        path_mct    = [path_mct,path_file];
                        path_mct    = path_mct(1:end-1);
                        obj.Rtomo   = fitsread(path_mct)';
                        obj.Rtomo   = obj.Rtomo(1:end-2,:);
                        %Wave front reconstruction
                        obj.waveFront = obj.Rtomo*reshape(obj.slopes(:,1:end-1,:),size(obj.Rtomo,2),[]);
                        % Get the tip-tilt
                        obj.tipTilt= zeros(2,nWfs,nExp);
                        for iWfs = 1:nWfs
                            obj.tipTilt(:,iWfs,:)= obj.Rtt*squeeze(obj.slopes(:,iWfs,:));
                        end
                    else
                        %Wave front reconstruction
                        obj.waveFront = obj.R*squeeze(obj.slopes(:,its,:));
                        obj.tipTilt = obj.Rtt*squeeze(obj.slopes(:,its,:));
                    end
                    
                    
                    %4\ Get Influence DM functions
                    bif      = gaussianInfluenceFunction(0.2);
                    nActu1D  = 8;
                    nRes     = 2*nActu1D+1;
                    dmSq     = deformableMirror(nActu1D,'modes',bif,'resolution',nRes);
                    x        = linspace(-1,1,8);
                    [x,y]    = meshgrid(x);
                    r        = hypot(x,y);
                    msk      = r<(1.0+1.4/7) & r>(0.285-1./7);
                    obj.dmIF = dmSq.modes.modes(:,msk);%unitless
                    
                    
                    %5\ Get the loop status
                    %5.1 Loop parameters
                    obj.holoop_lat  = 0.003;
                    obj.holoop_freq = cell2mat(val(strcmp(list,'FREQ')));
                    obj.holoop_gain = cell2mat(val(strcmp(list,'GAINDM')));
                    obj.ttloop_lat  = 0.003;
                    obj.ttloop_gain = cell2mat(val(strcmp(list,'GAINTX')));
                    obj.ttloop_freq = obj.holoop_freq;
                    
                    %5.2 Rejection transfer function
                    % ----------  Integrator tf
                    ho_num  = [-obj.holoop_gain 0];
                    ho_den  = [1 -1];
                    h_servo = tf(ho_num,ho_den,1/obj.holoop_freq);
                    
                    % Lag tf
                    delay  = obj.holoop_lat*obj.holoop_freq  + 1.05;
                    delta  = (delay - floor(delay));
                    h_lag  = tf('z',1/obj.holoop_freq)^(-1) + delta*tf('z',1/obj.holoop_freq)^(-2);
                    nF     = size(obj.waveFront,2);
                    floc           = logspace(-3,log10(0.5*obj.holoop_freq),nF/2);
                    obj.holoop_rtf = squeeze(bode(1/(1+h_servo*h_lag),2*pi*floc))';
                    
                    %Noise transfer function
                    obj.holoop_ntf = squeeze(bode(h_servo*h_lag/(1+h_servo*h_lag),2*pi*floc))';
                    obj.holoop_pn  = (trapz(floc,abs(obj.holoop_ntf).^2)*2/obj.holoop_freq);
                    
                    % --------------- TT Integrator tf
                    tt_num  = [-obj.ttloop_gain 0];
                    tt_den  = [1 -1];
                    ht_servo = tf(tt_num,tt_den,1/obj.ttloop_freq);
                    
                    % TT Lag tf
                    delay  = obj.ttloop_lat*obj.ttloop_freq;
                    delta  = (delay - floor(delay));
                    if delay > 1
                        ht_lag = tf('z',1/obj.ttloop_freq)^(-1) + delta*tf('z',1/obj.ttloop_freq)^(-2);
                    else
                        ht_lag = 1 + delta*tf('z',1/obj.ttloop_freq)^(-1);
                    end
                    nF = size(obj.tipTilt,2);
                    floc           = logspace(-3,log10(0.5*obj.ttloop_freq),nF/2);
                    obj.ttloop_rtf = squeeze(bode(1/(1+ht_servo*ht_lag),2*pi*floc))';
                    %Noise transfer function
                    obj.ttloop_ntf = squeeze(bode(ht_servo*ht_lag/(1+ht_servo*ht_lag),2*pi*floc))';
                    obj.ttloop_pn  = (trapz(floc,abs(obj.ttloop_ntf).^2)*2/obj.ttloop_freq);
                end
            end
        end
    end
    
end
