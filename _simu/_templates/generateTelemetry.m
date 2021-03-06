function trs = generateTelemetry(tel,atm,gs,sref,wfs,ts,nIter,varargin)
inputs = inputParser;
inputs.addRequired('tel', @(x) isa(x,'telescope'));
inputs.addRequired('atm', @(x) isa(x,'atmosphere'));
inputs.addRequired('gs',@(x) isa(x,'source'));
inputs.addRequired('sref',@(x) isa(x,'source'));
inputs.addRequired('wfs', @(x) isa(x,'shackHartmann'));
inputs.addRequired('ts', @(x) isa(x,'shackHartmann'));
inputs.addRequired('nIter', @isnumeric);
inputs.addParameter('gsRef', [],@(x) isa(x,'source'));
inputs.addParameter('S2Z', [],@isnumeric);
inputs.addParameter('P2Z', [],@isnumeric);
inputs.addParameter('frozenflow', false,@islogical);
inputs.addParameter('getZernike', false,@islogical);
inputs.addParameter('mmse', true,@islogical);
inputs.addParameter('ron', 0,@isnumeric);
inputs.addParameter('multicpu',true,@islogical);
inputs.parse(tel,atm,gs,sref,wfs,ts,nIter,varargin{:});

frozenflow  = inputs.Results.frozenflow;
mmse        = inputs.Results.mmse;
ron         = inputs.Results.ron;
S2Z         = inputs.Results.S2Z;
P2Z         = inputs.Results.P2Z;
multicpu    = inputs.Results.multicpu;
getZernike  = inputs.Results.getZernike;
gsRef       = inputs.Results.gsRef;

%% 1\ INIT
close all;
tel     = tools.duplicateTelescope(tel);
tel     = tel + atm; %attach atmosphere to telescope
sref    = sref.*tel*ts; %Define the light optical path for the Truth sensor;
gs      = gs.*tel*wfs; %Define the light optical path for the Truth sensor;

% Instantiate outputs
% Slopes time-vector
nSl     = size(wfs.slopes,1);
nGs     = wfs.lenslets.nArray;
tsSl    = zeros(nSl,nIter);
wfsSl   = zeros(nSl,nGs,nIter);
% WFSs camera
nPx     = ts.camera.resolution(1);
tsCam   = zeros(nPx,nPx,nIter);
wfsCam  = zeros(nPx,nPx*nGs,nIter);
% Optical phase difference map
opdTS   = zeros(tel.resolution,tel.resolution,nIter);
opdGS  = zeros(tel.resolution,tel.resolution,nGs,nIter);

%Enable readout and photon noise
if ron
    wfs.camera.readOutNoise = ron;
    wfs.camera.photonNoise = true;    
end

%% 2\ MULTI-CPU LOOP FOR GENERATING SIMULATED DATA SETS
addAttachedFiles(gcp,{'telescope.m','telescopeAbstract.m','shackHartmann.m','source.m'})


if multicpu && ~frozenflow
    % multi-cpu code
    for kIter=1:nIter
        %1\ Updating phase screens
        rngStream = RandStream('mt19937ar');
        draw(tel,rngStream);
       
        %2\ Propagating sources to the TS
        s2 = times(sref,tel);
        s2 = mtimes(s2,ts);
        %store TS slopes and pixels
        tsSl(:,kIter)       = ts.slopes;
        tsCam(:,:,kIter)    = ts.camera.frame;
        opdTS(:,:,kIter)    = s2.meanRmOpd;
        
        %3\ Propagating sources to the NGS WFS
        n2 = times(gs,tel);
        n2 = mtimes(n2,wfs);
        
        if ~isempty(gsRef)
            n3 = times(gsRef,tel);
            opdGS(:,:,:,kIter) = reshape([n3.meanRmOpd],tel.resolution,tel.resolution,nGs);
        end
        
        %store slopes and pixels
        wfsSl(:,:,kIter)    = wfs.slopes;
        wfsCam(:,:,kIter)   = wfs.camera.frame;
    end
    
else
    % single CPU code
    for kIter=1:nIter
        if ~mod(kIter,round(nIter/10))
            kIter
        end
        
        %1\ Updating phase screens
        if ~frozenflow
            draw(tel);
        else
            uplus(tel); % do not work in multi-cpu as the results a kIter depends on kIter-1
        end
        
        %2\ Propagating sources to the TS
        s2 = times(sref,tel);
        s2 = mtimes(s2,ts);
        %store TS slopes and pixels
        tsSl(:,kIter) = ts.slopes;
        tsCam(:,:,kIter) = ts.camera.frame;
        opdTS(:,:,kIter) = s2.meanRmOpd;
        
        %3\ Propagating sources to the NGS WFS
        n2 = times(gs,tel);
        n2 = mtimes(n2,wfs);
        %store slopes and pixels
        wfsSl(:,:,kIter) = wfs.slopes;
        wfsCam(:,:,kIter) = wfs.camera.frame;
        opdGS(:,:,:,kIter) = reshape([n2.meanRmOpd],tel.resolution,tel.resolution,nGs);
    end
end

trs = [];
trs.tsSl   = tsSl;
trs.tsCam  = tsCam;
trs.opdTS  = opdTS;
trs.wfsSl  = wfsSl;
trs.wfsCam = wfsCam;
trs.opdGS  = opdGS;

%% 3\ UPDATE RESULTS STRUCTURE
if ~isempty(S2Z)
    %1\ Wavefront errors
    trs.wfe = [];
    trs.wfe.uncorrected     = getWaveFrontErrorFromSlopes(tsSl,S2Z);
    trs.wfe.uncorrected_th  = sqrt(1.03*(tel.D/atm.r0)^(5/3))*atm.wavelength*1e9/2/pi;
    %2\ MMSE
    if mmse
        % note: with LGSs, we may have several sources at the same position
        % in the field to account for the sodium profile thickness; we take
        % the lowest ones
        [trs.Rmmse,trs.wfe.mmse_th] = getMMSE(tel,atm,gs(1:nGs),wfs,sref,S2Z);
        %trs.tomoSl                  = trs.Rmmse*reshape(trs.wfsSl,nSl*nGs,nIter);
        
        %iNoiseCovMat = wfs.theoreticalNoise(tel, atm, gs, sref,'computeInverseCovMat',true);%,'naParam',repmat([10000,90000],3,1),'lgsLaunchCoord',zeros(3,2));
        %trs.Rmmse    = sparseLMMSE(wfs,tel,atm,gs);%,'iNoiseVar',blkdiag(iNoiseCovMat{:}));
        trs.tomoSl   = trs.Rmmse*reshape(trs.wfsSl,nSl*nGs,nIter);
        trs.wfe.mmse = getWaveFrontErrorFromSlopes(trs.tsSl - trs.tomoSl,S2Z);
    end
    %3\ Zernike coefficients
    if getZernike && ~isempty(P2Z)
        % ON-AXIS
        trs.zer.ts_truth = P2Z * sref.waveNumber*reshape(opdTS,tel.resolution^2,nIter);
        trs.zer.ts_rec = S2Z * tsSl;
        
        % OFF-AXIS
        nZer = size(P2Z,1);
        trs.zer.wfs_truth = zeros(nGs,nZer,nIter);
        trs.zer.wfs_rec   = zeros(nGs,nZer,nIter);
        for n = 1:nGs
            trs.zer.wfs_truth(n,:,:) = P2Z * sref.waveNumber*reshape(squeeze(opdGS(:,:,n,:)),tel.resolution^2,nIter);
            trs.zer.wfs_rec(n,:,:)   = S2Z * squeeze(wfsSl(:,n,:));
        end
        % TOMOGRAPHY
        if mmse
           trs.zer.mmse = S2Z * trs.tomoSl; 
        end
    end
end


