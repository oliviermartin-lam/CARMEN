function [Cn2_d, alt_d, r0,date_bin, Cn2_upd,r0_upd,date_upd] = readScidarData(path,profileConfig,dt,nL_c,wvl,varargin)
inputs = inputParser;
inputs.addRequired('path', @ischar);
inputs.addRequired('profileConfig',  @ischar);
inputs.addRequired('dt', @isnumeric);
inputs.addRequired('nL_c',@isnumeric);
inputs.addRequired('wvl', @isnumeric);
inputs.addParameter('tupd', [], @isnumeric);
inputs.addParameter('display', false, @islogical);
inputs.parse(path,profileConfig,dt,nL_c,wvl,varargin{:});
tupd    = inputs.Results.tupd;
display = inputs.Results.display;


switch profileConfig
    case 'calm'
        file = '20140716.txt';
    case 'typical'
        file = '20141006.txt';
    case 'variable'
        file = '20151005.txt';
end

%For each night file there is a row per profile with format:
%YYYY-MM-DDTHH:MM:SS.ss,r0,seeing,coherenceTime,isoplanaticAngle,scintillationIndex,alt_0,cn2_0,windSpeed_0,windDirection_1,al_t1,cn2_1, windSpeed _1, windDirection _1….alt_n,cn2_n, windSpeed _n, windDirection _n

%% read
T     = readtable([path,file],'HeaderLines',0,'Delimiter',',');
nDate = numel(T{:,1});
alt   = 0:250:24750;
nL    = numel(alt);
Cn2   = zeros(nDate,nL);
loc_date = zeros(1,nDate);
for k=1:nDate
    tmp = cell2mat(T{k,1});
    loc_date(k) = str2double(tmp(12:13)) + str2double(tmp(15:16))/60 + str2double(tmp(18:19))/3600;
    if loc_date(k) > 12
        loc_date(k) = loc_date(k)-24;
    end
    for j=1:nL
        Cn2(k,j)  = T{k,8 + 4*(j-1)};
    end
end


%% Compressing the profile
nBin    = size(Cn2,1);
Cn2_c   = zeros(nBin,nL_c);
alt_d   = zeros(nBin,nL_c);

for k=1:nBin
    [Cn2_c(k,:),alt_d(k,:)] = eqLayers(Cn2(k,:),alt,nL_c);
    alt_d(k,alt_d(k,:)~=alt_d(k,:)) = 0;
end

%% Temporal decimation
dt_ori      = abs(median(diff(loc_date)))*60;
n_dec       = max(1,round(dt/dt_ori));
Cn2_d       = Cn2_c(1:n_dec:end,:);
date_bin    = loc_date(1:n_dec:end);

% r0 calculation
r0    = ((2*pi/wvl)^2 * 0.423*sum(Cn2_d,2)).^(-3/5);

if display
    figure;
    plot(date_bin,r0);
    ylabel('$r_0$ (m)','interpreter','latex','fontsize',18);
    xlabel('Time (h)','interpreter','latex','fontsize',18);
    set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
end

% Deriving the Cn2 profile at the MMSE reconstructo updates
if tupd
    n_upd    = max(1,round(tupd/dt_ori));
    Cn2_upd  = Cn2_c(1:n_upd:end,:);
    date_upd = loc_date(1:n_upd:end);
else
    Cn2_upd  = Cn2_c(1,:);
    date_upd = loc_date(1);
end
r0_upd   = ((2*pi/wvl)^2 * 0.423*sum(Cn2_upd,2)).^(-3/5);
