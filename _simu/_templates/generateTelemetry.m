function trs = generateTelemetry(tel,atm,ngs,sref,wfs,ts,nIter,varargin)
inputs = inputParser;
inputs.addRequired('tel', @(x) isa(x,'telescope'));
inputs.addRequired('atm', @(x) isa(x,'atmosphere'));
inputs.addRequired('ngs',@(x) isa(x,'source'));
inputs.addRequired('sref',@(x) isa(x,'source'));
inputs.addRequired('wfs', @(x) isa(x,'shackHartmann'));
inputs.addRequired('ts', @(x) isa(x,'shackHartmann'));
inputs.addRequired('nIter', @isnumeric);
inputs.addParameter('S2Z', [],@isnumeric);
inputs.addParameter('training', true,@islogical);
inputs.addParameter('mmse', true,@islogical);
inputs.addParameter('ron', 0,@isnumeric);
inputs.addParameter('multicpu',true,@islogical);
inputs.parse(tel,atm,ngs,sref,wfs,ts,nIter,varargin{:});

training    = inputs.Results.training;
mmse        = inputs.Results.mmse;
ron         = inputs.Results.ron;
S2Z         = inputs.Results.S2Z;
multicpu    = inputs.Results.multicpu;

%% 1\ INIT
close all;
tel     = tel + atm; %attach atmosphere to telescope
sref    = sref.*tel*ts; %Define the light optical path for the Truth sensor;
ngs     = ngs.*tel*wfs; %Define the light optical path for the Truth sensor;

% Instantiate outputs
% Slopes time-vector
nSl     = size(wfs.slopes,1);
nNgs    = numel(ngs);
tsSl    = zeros(nSl,nIter);
wfsSl   = zeros(nSl,nNgs,nIter);
% WFSs camera
nPx     = wfs.camera.resolution(1);
tsCam   = zeros(nPx,nPx,nIter);
wfsCam  = zeros(nPx,nPx*nNgs,nIter);
% Optical phase difference map
opdTS   = zeros(tel.resolution,tel.resolution,nIter);
opdNGS  = zeros(tel.resolution,tel.resolution,nNgs,nIter);

%Enable readout and photon noise
if ron
    wfs.camera.readOutNoise = ron;
    wfs.camera.photonNoise = true;    
end

%% 2\ MULTI-CPU LOOP FOR GENERATING SIMULATED DATA SETS
addAttachedFiles(gcp,{'telescope.m','telescopeAbstract.m','shackHartmann.m','source.m'})


if multicpu && training
    % multi-cpu code
    parfor kIter=1:nIter
        %1\ Updating phase screens
        draw(tel);
       
        %2\ Propagating sources to the TS
        s2 = times(sref,tel)
        s2 = mtimes(s2,ts)
        %store TS slopes and pixels
        tsSl(:,kIter) = ts.slopes;
        tsCam(:,:,kIter) = ts.camera.frame;
        opdTS(:,:,kIter) = s2.meanRmOpd;
        
        %3\ Propagating sources to the NGS WFS
        n2 = times(ngs,tel)
        n2 = mtimes(n2,wfs)
        %store slopes and pixels
        wfsSl(:,:,kIter) = wfs.slopes;
        wfsCam(:,:,kIter) = wfs.camera.frame;
        opdNGS(:,:,:,kIter) = reshape([n2.meanRmOpd],tel.resolution,tel.resolution,nNgs);
    end
    
else
    % single CPU code
    for kIter=1:nIter
        if ~mod(kIter,round(nIter/10))
            kIter
        end
        
        %1\ Updating phase screens
        if training
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
        n2 = times(ngs,tel);
        n2 = mtimes(n2,wfs);
        %store slopes and pixels
        wfsSl(:,:,kIter) = wfs.slopes;
        wfsCam(:,:,kIter) = wfs.camera.frame;
        opdNGS(:,:,:,kIter) = reshape([n2.meanRmOpd],tel.resolution,tel.resolution,nNgs);
    end
end

trs = [];
trs.tsSl   = tsSl;
trs.tsCam  = tsCam;
trs.opdTS  = opdTS;
trs.wfsSl  = wfsSl;
trs.wfsCam = wfsCam;
trs.opdNGS = opdNGS;

%% 3\ UPDATE RESULTS STRUCTURE
%1\ Wavefront errors
if ~isempty(S2Z)
    trs.wfe = [];
    trs.wfe.uncorrected     = getWaveFrontErrorFromSlopes_CANARY(tsSl);
    trs.wfe.uncorrected_th  = sqrt(1.03*(tel.D/atm.r0)^(5/3))*atm.wavelength*1e9/2/pi;
    %2\ MMSE
    if mmse
        [trs.Rmmse,trs.wfe.mmse_th] = getMMSE_CANARY(tel,atm,ngs,wfs,sref,S2Z);
        trs.tomoSl                  = trs.Rmmse*reshape(trs.wfsSl,nSl*nNgs,nIter);
        trs.wfe.mmse                = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl - trs.tomoSl);
    end
end


