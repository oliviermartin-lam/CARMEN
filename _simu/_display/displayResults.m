function displayResults(n,phLgs,phTs,wfsCam,tsCam,wfeLgs,wfeTs,wfeTomo,wfeRes,kIter,nIter,nExp,nPxWfs)

% display phase and pixels intensity
figure(n);
subplot(3,5,1)
imagesc(squeeze(phLgs(:,:,1,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
title('LGS 1','interpreter','Latex','fontSize',18);
subplot(3,5,2)
imagesc(squeeze(phLgs(:,:,2,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
title('LGS 2','interpreter','Latex','fontSize',18);
subplot(3,5,3)
imagesc(squeeze(phTs(:,:,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
title('TS','interpreter','Latex','fontSize',18);
subplot(3,5,4)
imagesc(squeeze(phLgs(:,:,3,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
title('LGS 3','interpreter','Latex','fontSize',18);
subplot(3,5,5)
imagesc(squeeze(phLgs(:,:,4,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
title('LGS 4','interpreter','Latex','fontSize',18);

subplot(3,5,6)
imagesc(squeeze(wfsCam(:,1:nPxWfs,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
subplot(3,5,7)
imagesc(squeeze(wfsCam(:,1+nPxWfs:2*nPxWfs,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
subplot(3,5,8)
imagesc(squeeze(tsCam(:,:,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
subplot(3,5,9)
imagesc(squeeze(wfsCam(:,1+2*nPxWfs:3*nPxWfs,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);
subplot(3,5,10)
imagesc(squeeze(wfsCam(:,1+3*nPxWfs:end,kIter)));
pbaspect([1,1,1])
set(gca,'XTick',[],'Ytick',[]);

subplot(3,5,11:12)
plot((1:kIter)*nExp/nIter,wfeLgs(1,1:kIter),'k');hold on;
plot((1:kIter)*nExp/nIter,wfeLgs(2,1:kIter),'b');
plot((1:kIter)*nExp/nIter,wfeLgs(3,1:kIter),'r');
plot((1:kIter)*nExp/nIter,wfeLgs(4,1:kIter),'m');
ylabel('LGS WFE (nm)','interpreter','Latex','fontSize',18);
xlabel('Ellapsed time (s)','interpreter','Latex','fontSize',18);
pbaspect([1.6,1,1])

subplot(3,5,14:15)
plot((1:kIter)*nExp/nIter,wfeTs(1:kIter),'k');hold on;
plot((1:kIter)*nExp/nIter,wfeTomo(1:kIter),'b');
plot((1:kIter)*nExp/nIter,wfeRes(1:kIter),'r--');
ylabel('TS WFE (nm)','interpreter','Latex','fontSize',18);
xlabel('Ellapsed time (s)','interpreter','Latex','fontSize',18);
legend({'Truth sensor','Tomography','Residual'},'interpreter','latex','FontSize',10,'Location','southwest');
pbaspect([1.6,1,1])

end