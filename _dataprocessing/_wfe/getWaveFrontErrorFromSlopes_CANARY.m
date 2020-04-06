function [wfe,bias,wfe_t] = getWaveFrontErrorFromSlopes_CANARY(slopes,varargin)
inputs = inputParser;
inputs.addRequired('slopes', @isnumeric);
inputs.addParameter('S2Z', [],@isnumeric);
inputs.parse(slopes,varargin{:});
S2Z = inputs.Results.S2Z;
%son must be provided as a matrix of size 72 x nIteration

% Manage the Zernike reconstruction matrix
if isempty(S2Z) 
    psTS = 0.219245;
    S2Z =  (0.5 * 4.2 * 1e9 * psTS / (3600.*180./pi)) * fitsread('GLOB_mrz_7x7.fits')'; %in nm
end

% Mean removal
smean = mean(slopes,2);
s_rmmean = bsxfun(@minus,slopes,smean);%subtract the temporal average

%Zernike decomposition
z_mean = S2Z*smean; % coefficient of the averaged centroids 35x1
z_rmmean = S2Z*s_rmmean; % mean-removed 35x2048

%Get wfe errors
wfe_t = sqrt(sum(z_rmmean.^2,1)); % sum of mean-removed Zernike coefficients
wfe = mean(wfe_t);
bias = sqrt(sum(z_mean.^2));% Mean value
