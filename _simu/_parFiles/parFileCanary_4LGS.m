%1\ TELESCOPE
D       = 4.2; %telescope diameter in meter
cobs 	= 0.285; % central obscuration ratio
Fe      = 150; % Phase screen upgrate rate in Hz


%2\ ATMOSPHERE 
%2\ ATMOSPHERE
atmName         = 'test1';

switch atmName
    
    case 'test1'
        photoAtm        = photometry.V0;      % photometric band to define atmospheric parameters V0->500nm
        windSpeed       = [7.5, 12.5, 15., 20.];       % wind speed vector in m/s
        windDirection	= [0, 11*pi/6, 3*pi/4, 4*pi/3];    % wind direction vector in rad
        L0           	= [25 25 25 25];         % Outer scale in meters
        r0           	= 0.16;                         % coherence lenght in meters at 500 nm
        fractionalR0 	= [65,15,10,10]/100.; % fractional weight of turbulent layers
        altitude        = [0, 4000., 10000., 15500.];
        
    case 'test2'
        %2\ ATMOSPHERE - Test 2
        photoAtm        = photometry.V0;      % photometric band to define atmospheric parameters V0->500nm
        windSpeed       = [7.5, 12.5, 15., 20.];       % wind speed vector in m/s
        windDirection	= [0, 11*pi/6, 3*pi/4, 4*pi/3];    % wind direction vector in rad
        L0           	= [25 25 25 25];         % Outer scale in meters
        r0           	= 0.12;                         % coherence lenght in meters at 500 nm
        fractionalR0 	= [45,15,30,10]/100.; % fractional weight of turbulent layers
        altitude        = [0, 2500., 4000., 13500.];
        
    case 'test3'
        %2\ ATMOSPHERE - Test 3
        photoAtm        = photometry.V0;      % photometric band to define atmospheric parameters V0->500nm
        windSpeed       = [10., 15., 17.5, 25.];       % wind speed vector in m/s
        windDirection	= [0, 11*pi/6, 3*pi/4, 4*pi/3];    % wind direction vector in rad
        L0           	= [25 25 25 25];         % Outer scale in meters
        r0           	= 0.085;                         % coherence lenght in meters at 500 nm
        fractionalR0 	= [80,5,10,5]/100.; % fractional weight of turbulent layers
        altitude        = [0, 6500., 10000., 15500.];
end

%3\ GSs
photoGs        = photometry.V;                % photometric band    550 nm
rGs            = [23, 23, 23, 23]*constants.arcsec2radian; % LGS asterism in phase B
dGs            = [0, pi/2, pi, 3*pi/2];
magGs          = [10, 10, 10, 10];
hGs            = 22.5e3; %altitude in meters
spotFwhm       = 1; % spot FWHM in arcsec
nLayerNa       = 2;
naWeight       = ones(1,nLayerNa); %sodium profile
naThickness    = 1e3;
viewPoint      = [0,0]; %launch from the secondary mirror

%4\ NGS WFS
% wfs geometry
nL              = 7; %\# 1D lenselts
nPx             = 16; %1D \# pixels per lenslet
nPxWfs          = nL*nPx;
d               = D/nL; %subaperture size
minLightRatio   = 0.5; % ratio of subap illumination to be valid
nyquistFlag     = false;
wfsPscale       = 0.47;% Pixel scale in arcsec for LGS WFS - Morris 2014
fovWfs          = wfsPscale*nPx; %WFS field of view in arcsec, must be close to 4"
% Calculate the number of pixels to get the closest pixel scale
lambdaOverd     = constants.radian2arcsec*photoGs.wavelength/d;
fovWfsinlod     = floor(fovWfs/lambdaOverd); % WFS fov in l/d units
Samp            = lambdaOverd/wfsPscale; % WFS sampling;
resTel          = nPx*nL*round(1/Samp); % \# pixels within the pupil
fovTel          = 2*max(rGs)*constants.radian2arcsec; % telescope fov ni arcsec.
ron             = 0.7; %Vidal et. al. 2014
wfsQE           = 0.1;% include detector quantum efficiency and throughput

%5\ LOOP
nIter           = 100;%number of simulated frames
nZern           = 35;

%6\ TRAINING
nScreens        = 5; %% TBC %%
hmin            = 10;%minimal altitude
hmax            = 15500;%maximal altitude
L0_training     = [25,25];
fractionalR0_training = [0.5,0.5];

%7\ OBS
profileConfig = 'typical'; %typical, calm or variable
nL_c = 6; % number of layers
dt   = 5; %in minutes