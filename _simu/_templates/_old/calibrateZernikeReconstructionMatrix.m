function S2Z = calibrateZernikeReconstructionMatrix(tel,wfs,sref,varargin)
inputs = inputParser;
inputs.addRequired('tel', @(x) isa(x,'telescope'));
inputs.addRequired('wfs', @(x) isa(x,'shackHartmann'));
inputs.addRequired('sref',@(x) isa(x,'source'));
inputs.addParameter('nModes', 200,@isnumeric);
inputs.addParameter('nThresholded', 1,@isnumeric);
inputs.parse(tel,wfs,sref,varargin{:});
nModes = inputs.Results.nModes;
nThresholded = inputs.Results.nThresholded;

%1\ Init
wfs.camera.readOutNoise = 0;
wfs.camera.photonNoise = false;
wfs.slopes = 0*wfs.slopes;

%2\ Define Zernike modes, piston-excluded
zer   = zernike(2:nModes+1,'resolution',tel.resolution);

%3\ Calibration
dm = deformableMirror(nModes,'modes',zer,'resolution',tel.resolution);
sref = sref.*tools.duplicateTelescope(tel); %duplicate the telescope so as to avoid including the atmosphere durinf the calibration process
dmCalib = calibration(dm,wfs,sref,sref.wavelength/40,1);
dmCalib.nThresholded = nThresholded;

% Zernike reconstruction matrix
S2Z = dmCalib.M;
