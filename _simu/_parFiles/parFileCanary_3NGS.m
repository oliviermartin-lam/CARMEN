\%1\ TELESCOPE
D 	= 4.2; %telescope diameter in meter
cobs 	= 0.285; % central obscuration ratio
Fe      = 150; % Phase screen upgrate rate in Hz


%2\ ATMOSPHERE
photoAtm	= photometry.V0;      % photometric band to define atmospheric parameters V0->500nm
windSpeed	= [7.5, 12.5, 15., 20.];       % wind speed vector in m/s
windDirection	= [0, 11*pi/6, 3*pi/4, 4*pi/3];    % wind direction vector in rad
L0           	= [25 25 25 25];         % Outer scale in meters
r0           	= 0.16;                         % coherence lenght in meters at 500 nm
fractionalR0 	= [65,15,10,10]/100.; % fractional weight of turbulent layers


%3\ NGSs
photoNgs= photometry.R;                % photometric band    640 nm
rNgs 	= [40.6 53 47.9]*constants.arcsec2radian;
dNgs 	= [pi/4 pi/2+pi/6,pi+pi/3];
magNgs  = [10.2 8.7 9.9];

%4\ NGS WFS
% wfs geometry
nL	= 7; %\# 1D lenselts
nPx     = 8; %1D \# pixels per lenslet
nPxWfs  = nL*nPx;
d       = D/nL; %subaperture size
minLightRatio  = 0.5; % ratio of subap illumination to be valid
wfsPscale = 0.21;% Pixel scale in arcsec
fovWfs    = wfsPscale*nPx; %WFS field of view in arcsec
% Calculate the number of pixels to get the closest pixel scale
lambdaOverd = constants.radian2arcsec*photoNgs.wavelength/d;
fovWfsinlod = floor(fovWfs/lambdaOverd); % WFS fov in l/d units
Samp        = lambdaOverd/wfsPscale; % WFS sampling;
resTel      = nPx*nL*round(1/Samp); % \# pixels within the pupil
fovTel      = 2*max(rNgs)*constants.radian2arcsec; % telescope fov ni arcsec.

%5\ LOOP
nExp = 1;%simulation time in seconds. %% TBC %%
nScreens=1; %% TBC %%
