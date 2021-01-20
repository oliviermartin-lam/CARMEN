%1\ TELESCOPE
D               = 39;       %telescope diameter in meter
cobs            = 0.28;    % central obscuration ratio
Fe              = 500;      % Phase screen upgrate rate in Hz

%2\ ATMOSPHERE - Test 1
photoAtm        = photometry.V0;                    % photometric band to define atmospheric parameters V0->500nm
windSpeed       = [7.5, 12.5, 15., 20.];            % wind speed vector in m/s
windDirection	= [0, 11*pi/6, 3*pi/4, 4*pi/3];     % wind direction vector in rad
L0           	= [25 25 25 25];                    % Outer scale in meters
r0           	= 0.16;                             % coherence lenght in meters at 500 nm
fractionalR0 	= [65,15,10,10]/100.;               % fractional weight of turbulent layers
altitude        = [0, 4000., 10000., 15500.];


%3\ NGSs
photoNgs        = photometry.Na;                         % photometric band    640 nm
rNgs            = [70.0, 70.0, 70.0, 70.0, 70.0, 70.0]*constants.arcsec2radian; % A47 asterism in the CANARY nomenclature
dNgs            = [0 , 60, 120, 180, 240, 300];
magNgs          = [10, 10, 10, 10, 10, 10];

%4\ NGS WFS
% wfs geometry
nL              = 74;            %\# 1D lenselts
nPx             = 6;           %1D \# pixels per lenslet
nPxWfs          = nL*nPx;
d               = D/nL;         %subaperture size
minLightRatio   = 0.5;          % ratio of subap illumination to be valid
wfsPscale       = 0;         % Pixel scale in arcsec

if wfsPscale == 0
    resTel          = nPx*nL/2;
else
    fovWfs          = wfsPscale*nPx; %WFS field of view in arcsec
    
    % Calculate the number of pixels to get the closest pixel scale
    lambdaOverd     = constants.radian2arcsec*photoNgs.wavelength/d;
    fovWfsinlod     = floor(fovWfs/lambdaOverd);            % WFS fov in l/d units
    Samp            = lambdaOverd/wfsPscale;                % WFS sampling;
    resTel          = nPx*nL*round(1/Samp);                 % \# pixels within the pupil
end

fovTel          = 2*max(rNgs)*constants.radian2arcsec;  % telescope fov ni arcsec.
ron             = 0.2; 	
wfsQE           = 0.5;                                  % include detector quantum efficiency and throughput
%5\ LOOP
nIter           = 100;                                  %number of simulated frames

%6\ TRAINING
nScreens        = 5; %% TBC %%
hmin            = 10;%minimal altitude
hmax            = 20000;%maximal altitude
L0_training     = [25,25];
fractionalR0_training = [0.5,0.5];
