ticID = tic;

% add a postfix description to the folder name
if isfield( CONFIG.kernels, 'namePostfix' ) 
	RESULTS_path = [ 'Results_' CONFIG.kernels.namePostfix ];
else
	RESULTS_path = 'Results';
end

[~,tracking,~] = fileparts(CONFIG.TRACKING_path);
fprintf( '\n-> saving results to "%s/%s/Tracking/%s/%s/":\n', CONFIG.protocol, CONFIG.subject, tracking, RESULTS_path );
clear tracking

RESULTS_path = fullfile( CONFIG.TRACKING_path, RESULTS_path );
[~,~,~] = mkdir( RESULTS_path ); 

%% Load data from files for cross-validation:
nii1 = load_untouch_nii(fullfile(CONFIG.DATA_path, 'DWI.nii'));
nii2 = load_untouch_nii(fullfile(CONFIG.DATA_path, 'life_demo_scan2_subject1_b2000_150dirs_stanford.nii'));
bvals1 = load(fullfile(CONFIG.DATA_path, 'DWI.bvals'));
bvals2 = load(fullfile(CONFIG.DATA_path, 'life_demo_scan2_subject1_b2000_150dirs_stanford.bvals'));

%% Calculate S0, diffusion-weighted and relative signals:
s01 = mean(nii1.img(:, :, :, bvals1==0), 4);
s02 = mean(nii2.img(:, :, :, bvals2==0), 4);
diffusion_weighted1 = nii1.img(:, :, :, bvals1~=0);
diffusion_weighted2 = nii2.img(:, :, :, bvals2~=0);
relative_signal1 = double(diffusion_weighted1)./repmat(s01, 1, 1, 1,...
                    size(diffusion_weighted1, 4));
relative_signal2 = double(diffusion_weighted2)./repmat(s02, 1, 1, 1,...
                    size(diffusion_weighted2, 4));

%%
% configuration and results
% -------------------------
fprintf( '\t- configuration\n');
save( fullfile(RESULTS_path,'CONFIG.mat'), 'CONFIG' )
print( gcf, fullfile( RESULTS_path, sprintf('summary.png') ), '-dpng' );


% map of voxelwise errors
% -----------------------
fprintf( '\t- voxelwise map of fitting errors\n');

% CROSS-VALIDATION! Use the testing signal from the second measurement
% to calculate the prediction errors:
y_mea = double(diffusion_weighted2);

y_est = zeros( niiSIGNAL.hdr.dime.dim(2:5) );
y_est( DICTIONARY.MASKidx ) = A * CONFIG.OPTIMIZATION.x;
if ( CONFIG.useReference )
	y_mea = y_mea(2:end,:,:,:);
	y_est = y_est(2:end,:,:,:);
end

%%

% Permute this for comparison with the data, as saved in the files (bvals
% dimension last):
y_est = permute(y_est, [2, 3, 4, 1]);

% Add back the mean and multiply by the non diffusion-weighted signal
% As in: https://github.com/nipy/dipy/blob/master/dipy/tracking/life.py#L564
if CONFIG.doDemean
    mean_relative_signal1 = mean(relative_signal1, 4);
    y_est = y_est + repmat(mean_relative_signal1, 1, 1, 1, size(y_est, 4));
end

% Multiply by S0, to compare in scanner units:
y_est = y_est .* repmat(s01, 1, 1, 1, size(y_est, 4));


%%
% Reconstruction NRMSE
niiERR = niiSIGNAL;
niiERR.img = squeeze( sqrt(sum((y_mea-y_est).^2,4) ./ sum(y_mea.^2,4)) );
niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
niiERR.hdr.dime.datatype = 16;
niiERR.hdr.dime.bitpix   = 32;
niiERR.hdr.dime.dim(1:5) = [ 3 niiERR.hdr.dime.dim(3:5) 1 ];
niiERR.hdr.dime.pixdim(2:5) = niiERR.hdr.dime.pixdim(3:6);
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = 1;
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = 1;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_NRMSE.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('NRMSE'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_NRMSE.png') ), '-dpng' );
close( figure(99) )

%%
% Reconstruction RMSE
niiERR.img = squeeze( sqrt(mean((y_mea-y_est).^2,4)) );
niiERR.img( DICTIONARY.MASK==0 ) = 0;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = max( niiERR.img(:) );
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_RMSE.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('RMSE'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_RMSE.png') ), '-dpng' );
close( figure(99) )

%%
% Reconstruction SNR
niiERR.img = squeeze( 10*log10( mean(y_mea.^2,1) ./ mean((y_mea-y_est).^2,1) ) );
%niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = max( niiERR.img(:) );
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_SNR.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('SNR [dB]'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_SNR.png') ), '-dpng' );
close( figure(99) )

% Volume fractions
% ----------------
fprintf( '\t- volume fractions\n' );

if ( CONFIG.kernels.doNormalize )
	xW = CONFIG.OPTIMIZATION.x ./ [ reshape(repmat(KERNELS.wmr_norm,A.nF,1),1,[]) reshape(repmat(KERNELS.wmh_norm,A.nE,1),1,[]) reshape(repmat(KERNELS.iso_norm,A.nV,1),1,[]) ]';
else
	xW = CONFIG.OPTIMIZATION.x;
end

niiERR.hdr.dime.dim(1) = 4;
niiERR.hdr.dime.dim(5) = 3;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = 1;
niiERR.hdr.dime.cal_min = niiERR.hdr.dime.glmin;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
niiERR.img = zeros(niiERR.hdr.dime.dim(2:5),'single');

% isotropic
fprintf( '\t\t- isotropic\n' );
tmp  = zeros( DICTIONARY.nV, 1 );
base = A.nF*A.nR + A.nE*A.nT;
for s = 1:numel(DICTIONARY.ISO.v)
	v = DICTIONARY.ISO.v(s)/KERNELS.nS + 1;
	xx = xW(base+s : A.nV : base+A.nV*A.nI);
	for k = 1:A.nI
		tmp(v) = tmp(v) + xx(k);
	end
end
IMG  = zeros( DICTIONARY.dim );
IMG( DICTIONARY.MASK>0 ) = tmp;
niiERR.img( :,:,:,3 ) = IMG;

% extra-axonal
if numel(KERNELS.wmh) > 0
	fprintf( '\t\t- extra-axonal\n' );
    tmp  = zeros( DICTIONARY.nV, 1 );
    base = A.nF*A.nR;
	for s = 1:numel(DICTIONARY.EC.v)
		v = DICTIONARY.EC.v(s)/KERNELS.nS + 1;
		xx = xW( base+s: A.nE : base+A.nE*A.nT);
		for k = 1:A.nT
			tmp(v) = tmp(v) + xx(k);
		end
    end
    IMG  = zeros( DICTIONARY.dim );
    IMG( DICTIONARY.MASK>0 ) = tmp;
	niiERR.img( :,:,:,2 ) = IMG;
end

% intra-axonal
fprintf( '\t\t- intra-axonal\n' );
tmp   = zeros( DICTIONARY.nV, 1 );
tmp2  = zeros( DICTIONARY.nV, 1 );
for s = 1:DICTIONARY.IC.n
	f = DICTIONARY.IC.fiber(s) + 1;
	v = DICTIONARY.IC.v(s)/KERNELS.nS + 1;
	for k = 1:A.nR
		w = xW(f + (k-1)*A.nF) * DICTIONARY.IC.len(s);
		tmp(v) = tmp(v) + w;
	end
	tmp2(v) = tmp2(v) + DICTIONARY.IC.len(s);
end

IMG  = zeros( DICTIONARY.dim );
IMG( DICTIONARY.MASK>0 ) = tmp;
niiERR.img( :,:,:,1 ) = IMG;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'vf.nii') )

IMG  = zeros( DICTIONARY.dim );
IMG( DICTIONARY.MASK>0 ) = tmp2;
niiERR.img( :,:,:,1 ) = IMG;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'vf__raw.nii') )


fprintf( '   [ %.2f seconds ]\n', toc(ticID) );

clear ticID RESULTS_path niiERR y_est y_mea tmp tmp2 IMG w xW xx f v s k base

