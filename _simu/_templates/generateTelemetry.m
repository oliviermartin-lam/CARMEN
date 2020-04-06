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
inputs.parse(tel,atm,ngs,sref,wfs,ts,nIter,varargin{:});

training = inputs.Results.training;
mmse    = inputs.Results.mmse;
ron = inputs.Results.ron;
S2Z = inputs.Results.S2Z;

%% 1\ INIT
close all;
tel = tel + atm; %attach atmosphere to telescope
sref = sref.*tel*ts; %Define the light optical path for the Truth sensor;
ngs = ngs.*tel*wfs; %Define the light optical path for the Truth sensor;

% Instantiate outputs
trs = [];
% Slopes time-vector
nSl = size(wfs.slopes,1);
nNgs = numel(ngs);
trs.tsSl = zeros(nSl,nIter);
trs.wfsSl = zeros(nSl,nNgs,nIter);
% WFSs camera
nPx = wfs.camera.resolution(1);
trs.tsCam = zeros(nPx,nPx,nIter);
trs.wfsCam = zeros(nPx,nPx*nNgs,nIter);
% Optical phase difference map
trs.opdTS = zeros(tel.resolution,tel.resolution,nIter);
trs.opdNGS = zeros(tel.resolution,tel.resolution,nNgs,nIter);

%Enable readout and photon noise
if ron
    wfs.camera.readOutNoise = ron;
    wfs.camera.photonNoise = true;    
end

%% 2\ LOOP
for kIter=1:nIter
    if ~mod(kIter,round(nIter/10))
        kIter
    end
    
    %1\ Updating phase screens
    if training
        draw(tel);
    else
        +tel;
    end
        
    %2\ Propagating sources to the TS
    +sref;
    %store TS slopes and pixels
    trs.tsSl(:,kIter) = ts.slopes;
    trs.tsCam(:,:,kIter) = ts.camera.frame;
    trs.opdTS(:,:,kIter) = sref.meanRmOpd;
     
    %3\ Propagating sources to the NGS WFS
    +ngs;
    %store slopes and pixels
    trs.wfsSl(:,:,kIter) = wfs.slopes;
    trs.wfsCam(:,:,kIter) = wfs.camera.frame;
    for j=1:nNgs
        trs.opdNGS(:,:,j,kIter) = ngs(j).meanRmOpd;
    end
end

%% 3\ UPDATE RESULTS STRUCTURE
%1\ Wavefront errors
trs.wfe = [];
trs.wfe.uncorrected = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl);
trs.wfe.uncorrected_th = sqrt(1.03*(tel.D/atm.r0)^(5/3))*atm.wavelength*1e9/2/pi;
%2\ MMSE
if mmse
    [trs.Rmmse,trs.wfe.mmse_th] = getMMSE_CANARY(tel,atm,ngs,wfs,sref,S2Z);
    trs.tomoSl = trs.Rmmse*reshape(trs.wfsSl,nSl*nNgs,nIter);
    trs.wfe.mmse = getWaveFrontErrorFromSlopes_CANARY(trs.tsSl - trs.tomoSl);     
end


