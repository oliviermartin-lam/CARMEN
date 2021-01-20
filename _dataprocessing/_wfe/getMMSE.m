function [Rmmse,wfeTomo,Cee] = getMMSE(tel,atm,ngs,wfs,sref,S2Z)

if ~any([atm.layer.altitude])
    cond = 1e3;
else
    cond = 30;
end


covModel = slopesCovarianceModel(tel,atm,wfs,[sref,ngs],'Projector',S2Z,'isFFT',false);
ps = wfs.camera.pixelScale*constants.radian2arcsec;
covMat = covModel.getCovarianceMatrix()/ps^2; %in pixel^2
Rmmse = covModel.getMMMSEreconstructor(covMat,cond); % MMSE reconstructor
Cee  = covModel.getErrorCovarianceMatrix(covMat,Rmmse);
wfeTomo =  covModel.getWaveFrontError(Cee); %Cee in pixel
