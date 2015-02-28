% Load the WM mask: 
niiMASK = load_untouch_nii(fullfile(CONFIG.TRACKING_path,'dictionary_mask.nii'));

%% Compare errors between LiFE and COMMIT 

niiERR_L = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_NRMSE.nii') );
niiERR_C = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_NRMSE.nii') );

% plot the NRMSE of the COMMIT model
figure(2), clf, imagesc( rot90(squeeze(100*niiERR_C.img(:,70,:))), [0 100] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = 100*niiERR_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.1f%% +/- %.1f%%', mean(yL), std(yL) ))

% direct comparison of the NRMSE of LiFE and COMMIT
figure(3), clf, hold on
x = linspace(0,100,60);
yL = hist( 100*niiERR_L.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
yC = hist( 100*niiERR_C.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on, axis tight
xlabel( 'NRMSE [%]' ), ylabel( 'percentage of voxels' )
pbaspect([1 1 1])
legend('LiFE','COMMIT')
title('Error distributions')

% voxelwise comparison of the NRMSE of LiFE and COMMIT
figure(4), clf, hold on
yL = 100*niiERR_L.img( niiMASK.img>0 );
yC = 100*niiERR_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 100], [0 100], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 100 0 100])
xlabel( 'NRMSE [%] with LiFE' ), ylabel( 'NRMSE [%] with COMMIT' )
title('Error scatterplot')

% plot the NRMSE of the LiFE model
figure(1), clf, imagesc( rot90(squeeze(100*niiERR_L.img(:,70,:))), [0 100] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = 100*niiERR_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.1f%% +/- %.1f%%', mean(yL), std(yL) ))


niiERR_L  = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_RMSE.nii') );
niiERR_C  = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_RMSE.nii') );

% plot the RMSE of the LiFE model
figure(5), clf, imagesc( rot90(squeeze(niiERR_L.img(:,70,:))), [0 200] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = niiERR_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.1f +/- %.1f', nanmean(yL), nanstd(yL) ))
saveas(gcf,'RESULTS_Fig5.png')

% plot the RMSE of the COMMIT model
figure(6), clf, imagesc( rot90(squeeze(niiERR_C.img(:,70,:))), [0 200] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = niiERR_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.1f +/- %.1f', mean(yL), std(yL) ))
saveas(gcf,'RESULTS_Fig6.png')

% direct comparison of the RMSE of LiFE and COMMIT
figure(7), clf, hold on
x = linspace(0,300,100);
yL = hist( niiERR_L.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
yC = hist( niiERR_C.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on, axis tight
xlabel( 'RMSE [raw signal units]' ), ylabel( 'percentage of voxels' )
pbaspect([1 1 1])
legend('LiFE','COMMIT')
title('Error distributions')
saveas(gcf,'RESULTS_Fig7.png')

% voxelwise comparison of the RMSE of LiFE and COMMIT
figure(8), clf, hold on
yL = niiERR_L.img( niiMASK.img>0 );
yC = niiERR_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 260], [0 260], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 260 0 260])
xlabel( 'RMSE [raw signal units] with LiFE' ), ylabel( 'RMSE [raw signal units] with COMMIT' )
title('Error scatterplot')
saveas(gcf,'RESULTS_Fig8.png')