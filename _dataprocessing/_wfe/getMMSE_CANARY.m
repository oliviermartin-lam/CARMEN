function Rmmse = getMMSE_CANARY(tel,atm,ngs,wfs,sref)

covModel = slopesCovarianceModel(tel,atm,wfs,[sref,ngs]);
covMat = covModel.getCovarianceMatrix();
Rmmse = covModel.getMMMSEreconstructor(covMat); % MMSE reconstructor NGS+LGS slopes -> TS slopes

