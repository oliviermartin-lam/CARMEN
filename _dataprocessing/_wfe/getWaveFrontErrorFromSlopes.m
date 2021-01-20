function [wfe,bias,wfe_t] = getWaveFrontErrorFromSlopes(slopes,S2Z)
inputs = inputParser;
inputs.addRequired('slopes', @isnumeric);
inputs.addRequired('S2Z',@isnumeric);
inputs.parse(slopes,S2Z);
% Mean removal
smean       = mean(slopes,2);
s_rmmean    = bsxfun(@minus,slopes,smean);%subtract the temporal average
%Zernike decomposition
z_mean      = S2Z*squeeze(smean); % coefficient of the averaged centroids 35x1
z_rmmean    = S2Z*s_rmmean; % mean-removed 35x2048
%Get wfe errors
wfe_t       = sqrt(sum(z_rmmean.^2,1)); % sum of mean-removed Zernike coefficients
wfe         = mean(wfe_t);
bias        = sqrt(sum(z_mean.^2));% Mean value
