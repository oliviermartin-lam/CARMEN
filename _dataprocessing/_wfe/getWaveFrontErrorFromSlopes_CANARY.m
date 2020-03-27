function [wfe,bias] = getWaveFrontErrorFromSlopes_CANARY(sdiff)

%son must be provided as a matrix of size 72 x nIteration
%1\ Get the reconstructor
%volt2adu  = 2^15/10.;
%volt2nm= 0.78*1e3;
%adu2volt  = 1./volt2adu;
%tonm = volt2nm*adu2volt;
%R = fitsread('mc_2013-07-22_17h31m18s_TS.fits')';



R = fitsread('zmi_7x7.fits');
smean = mean(sdiff,2);
s_rmmean = bsxfun(@minus,sdiff,smean);%subtract the temporal average

z_mean = R*smean; % coefficient of the averaged centroids 35x1
z_rmmean = R*s_rmmean; % mean-removed 35x2048

wfe = sqrt(sum(z_rmmean.^2,1)); % sum of mean-removed Zernike coefficients
bias = sqrt(sum(z_mean.^2));% Mean value
