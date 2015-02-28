%%
addpath(genpath('~/source/COMMIT/'))

%%
COMMIT_Setup

%%
CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean = true;
DATA_Load

CONFIG.kernels.namePostfix  = 'LIFE';
CONFIG.kernels.d            = 1.7;
CONFIG.kernels.Rs           = [ 0 ];
CONFIG.kernels.ICVFs        = [  ];
CONFIG.kernels.dISOs        = [  ];

KERNELS_CreateFolderForHighResolutionKernels( CONFIG );
KERNELS_PrecomputeRotationMatrices();
KERNELS_StickZeppelinBall_Generate( CONFIG );
KERNELS_ActiveAx_RotateAndSave( CONFIG );

KERNELS = KERNELS_Load( CONFIG );
KERNELS_ProcessAtoms

%%
CONFIG.TRACKING_path        = fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup

%% 
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults