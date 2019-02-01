function [model_output, diff, R]= getModelDataDiff_realLGM_SST_Anomalies_d18OAnomalies(modelPath,ensSize,targetPathSST,targetPathD18O)
% Written by C.Breitkreuz (last modified 31.01.2019)
% getModelDataDiff_realLGM_SST_Anomalies_d18OAnomalies reads model output,
% computes simulated d18Oc and computes model-data difference.

% Input:
% modelPath - path to model output
% ensSize - total number of ensemble members 
% targetPathSST - path to SST anomalies data
% targetPathD18O - path to d18Oc anomalies data

% Output:
% model_output 
% diff - model-data misfif
% R - observation error covariance matrix


numObsFields = 3; % annual SST, planktonic and benthic foraminifera

PRef = 0;

%% Read model grid - LGM

gridDirectory = '/m/wrk3/cbreitkreuz/MITgcm_exp/global_ocean.cs32x15_LGM/run_LGM_wiso_2/mnc_output_all';


grids = rdmnc(fullfile(gridDirectory, 'grid.*'));

hFacC = grids.HFacC;
rC_LGM= abs(grids.RC);

% Determine array sizes
Nx = size(hFacC, 1);
Ny = size(hFacC, 2);
Nr = size(hFacC, 3);

mask_LGM = hFacC;
mask_LGM(mask_LGM > 0) =  1;
mask_LGM(mask_LGM==0) = nan;

%% Read Late holocene grid (C.Breitkreuz et al., 2018)

gridDirectory = '/m/wrk3/cbreitkreuz/MITgcm_exp/global_ocean.cs32x15_newMethod/';

grids = rdmnc(fullfile(gridDirectory, 'grid.*'));
hFacC = grids.HFacC;
rC_LH= abs(grids.RC);

mask_LH = hFacC;
mask_LH(mask_LH > 0) =  1;
mask_LH(mask_LH==0) = nan;

mask = mask_LH .* mask_LGM;

%% Read proxy data - SST LGM-L anomalies (on model grid)


targetFileSST_Annual = 'MARGO_SST_Anomalies_Annual_cs.bin';
targetFileSigma_Annual = 'MARGO_SST_Anomalies_SigmaAnnual_cs.bin';

SST_Annual_Anomalies = mit_readfield(fullfile(targetPathSST,targetFileSST_Annual), [192 32 15], 'real*8') .* mask;
SIGMA_Annual = mit_readfield(fullfile(targetPathSST,targetFileSigma_Annual), [192 32 15], 'real*8') .* mask;


SST_Annual_Anomalies = SST_Annual_Anomalies(:,:,1);
SIGMA_Annual = SIGMA_Annual(:,:,1);

SST_Annual_Anomalies(SIGMA_Annual==0)=nan;
SIGMA_Annual(SIGMA_Annual==0)=nan;


%% Read proxy data - d18O LGM-LH anomalies

targetFile_Planktonicd18O = 'LGM_PLANKTONIC_D18OAnom_DATA_cs.bin';
targetFileSigma_Planktonicd18O = 'LGM_PLANKTONIC_D18OAnom_ONESIGMA_cs.bin';

targetFile_Benthicd18O = 'LGM_BENTHIC_D18OAnom_cs.bin';
targetFileSigma_Benthicd18O = 'LGM_BENTHIC_D18OAnom_ONESIGMA_cs.bin';


d18O_calciteAnomalie_Planktonic = mit_readfield(fullfile(targetPathD18O,targetFile_Planktonicd18O ), [192 32 15 ], 'real*8') .* mask;
SIGMA_calciteAnomalie_Planktonic = mit_readfield(fullfile(targetPathD18O,targetFileSigma_Planktonicd18O), [192 32 15], 'real*8') .* mask;

d18O_calciteAnomalie_Planktonic = d18O_calciteAnomalie_Planktonic(:,:,1);
SIGMA_calciteAnomalie_Planktonic = SIGMA_calciteAnomalie_Planktonic(:,:,1);

d18O_calciteAnomalie_Benthic = mit_readfield(fullfile(targetPathD18O,targetFile_Benthicd18O ), [192 32 15], 'real*8') .* mask;
SIGMA_calciteAnomalie_Benthic = mit_readfield(fullfile(targetPathD18O,targetFileSigma_Benthicd18O), [192 32 15], 'real*8') .* mask;


d18O_calciteAnomalie_Planktonic(SIGMA_calciteAnomalie_Planktonic==0)=nan;
SIGMA_calciteAnomalie_Planktonic(SIGMA_calciteAnomalie_Planktonic==0)=nan;

d18O_calciteAnomalie_Benthic(SIGMA_calciteAnomalie_Benthic==0)=nan;
SIGMA_calciteAnomalie_Benthic(SIGMA_calciteAnomalie_Benthic==0)=nan;


%% Read d18O Late Holocene data set (C.Breitkreuz et al.,2018)

modelDirectory_LH ='/m/wrk3/cbreitkreuz/MITgcm_exp/ad_global_ocean.cs32x15/MOD_newCost/run_ad_wiso_400_1_iter13_forwardTest_400years_2/mnc_output_all';

% Read d18O Late Holocene
vIter = rdmnc(fullfile(modelDirectory_LH, 'tracerDiag.*'), 'iter');

iteration =  vIter.iter(end);
inputField_ptracer     = rdmnc(fullfile(modelDirectory_LH, 'tracerDiag.*'),iteration);

H2O16 = inputField_ptracer.TRAC01 .* mask;
H2O18 = inputField_ptracer.TRAC02 .* mask;

% Calculate d180 and correct from VSMOV to VPDB (- 0.27)
model_d18Osw_LH = (((H2O18./H2O16)./(1/498.7))-1) * 1000 -0.27 ;


% Read according temperature field
vIter = rdmnc(fullfile(modelDirectory_LH, 'dynDiag.*'), 'iter');


iteration =  vIter.iter(end);
inputField    = rdmnc(fullfile(modelDirectory_LH, 'dynDiag.*'),iteration,'THETA','SALT');

model_theta_LH = inputField.THETA .* mask;
model_salt_LH = inputField.SALT .* mask;


% Compute in-situ tempereature from potential temp. (THETA)
model_temp_LH = nan(Nx,Ny,Nr);
for k = 1:Nr
    model_temp_LH(:,:,k) = potential2insituT(model_salt_LH(:,:,k),model_theta_LH(:,:,k), PRef, rC_LH(k));
end


% Compute d18O calcite from it, according to Shakleton 
model_d18O_LH_calcite = model_d18Osw_LH + 21.9 - sqrt(310.61 + 10 .* model_temp_LH);


%% Start loop over ensemble members


model_theta_LGMAnnual = nan(Nx,Ny,Nr,ensSize);
model_salt_LGMAnnual = nan(Nx,Ny,Nr,ensSize);
model_SST_Annual_Anomalies = nan(Nx,Ny,ensSize);

model_d18Osw_LGM = zeros(Nx,Ny,Nr,ensSize);
model_d18O_swAnomalie= zeros(Nx,Ny,Nr,ensSize);

model_d18O_LGM_calcite= zeros(Nx,Ny,Nr,ensSize);
model_d18O_calciteAnomalie= zeros(Nx,Ny,Nr,ensSize);

% Read vector of iterations and find last time slice

if ensSize==1
    vIterYearly  = rdmnc(fullfile(fullfile(modelPath, '/mnc_output_all'), 'dynDiag.*'), 'iter');
else
    vIterYearly  = rdmnc(fullfile(fullfile(modelPath, '/mnc_output_all_1'), 'dynDiag.*'), 'iter');
end

itYearly = vIterYearly.iter(end);

diff = nan(Nx,Ny,Nr,numObsFields,ensSize);
model_output = nan(Nx,Ny,Nr,numObsFields,ensSize);

 model_temp_LGM  = nan(Nx,Ny,Nr,ensSize);
 
for ensembleMember = 1:ensSize
    
    fprintf(['Reading ensemble member ',num2str(ensembleMember),'\n'])
    
    %% Set directory name according to ensemble member/particle
    if ensSize==1
        modelExperiment = '/mnc_output_all';
    else
        modelExperiment = ['/mnc_output_all_',num2str(ensembleMember)];
    end
    
    modelDirectory  = fullfile(modelPath, modelExperiment);
    
    %% Read model temperature

    inputFieldAnnual = rdmnc(fullfile(modelDirectory, 'dynDiag.*'),'THETA','SALT', itYearly);
    
    model_theta_LGM_temp= inputFieldAnnual.THETA;
    model_salt_LGM_temp= inputFieldAnnual.SALT;
    
    model_theta_LGMAnnual(:,:,:,ensembleMember) = model_theta_LGM_temp .* mask;
    model_salt_LGMAnnual(:,:,:,ensembleMember) = model_salt_LGM_temp .* mask;

    
    % Compute in-situ tempereature from potential temp. (THETA)
    for k = 1:Nr
        model_temp_LGM (:,:,k,ensembleMember) = potential2insituT(model_salt_LGMAnnual(:,:,k,ensembleMember),model_theta_LGMAnnual(:,:,k,ensembleMember), PRef, rC_LGM(k));
    end

    
    %% Read model d18O

    % Read data structure from file
    inputField_ptracer     = rdmnc(fullfile(modelDirectory, 'tracerDiag.*'),'TRAC01','TRAC02',itYearly);
    
    H2O16_LGM = inputField_ptracer.TRAC01 .* mask;
    H2O18_LGM = inputField_ptracer.TRAC02 .* mask;
    
    % Calculate d180 and correct from VSMOV to VPDB (-0.27)and + 1.1 permil
    % correction for LGM ice sheet contribution
    model_d18Osw_LGM(:,:,:,ensembleMember) = (((H2O18_LGM./H2O16_LGM)./(1/498.7))-1) * 1000 - 0.27 + 1.1;
    
    % Compute d18O calcite, according to Shakleton
    model_d18O_LGM_calcite(:,:,:,ensembleMember)  = model_d18Osw_LGM(:,:,:,ensembleMember)  + 21.9 - sqrt(310.61 + 10 .* model_temp_LGM(:,:,:,ensembleMember) );

    
    %% Compute LGM-LH anomalies

    model_d18O_swAnomalie(:,:,:,ensembleMember)  = model_d18Osw_LGM(:,:,:,ensembleMember)  - model_d18Osw_LH;
    model_d18O_calciteAnomalie(:,:,:,ensembleMember)  = model_d18O_LGM_calcite(:,:,:,ensembleMember)  - model_d18O_LH_calcite;
    
    model_SST_Annual_Anomalies(:,:,ensembleMember) = model_temp_LGM(:,:,1,ensembleMember) - model_temp_LH(:,:,1) ;
    

    %% Model-output
    
    % Field 1, SST Anomalies
    field = 1;
    model_output(:,:,1,field,ensembleMember) = model_SST_Annual_Anomalies(:,:,ensembleMember); % surface level
    
    % Field 2, d18O calcite anomalie, planktonic
    field = 2;
    model_output(:,:,1,field,ensembleMember) = model_d18O_calciteAnomalie(:,:,1,ensembleMember); % 192 x 32
    
     % Field 3, d18O calcite anomalie, benthic
    field = 3;
    model_output(:,:,2:15,field,ensembleMember) = model_d18O_calciteAnomalie(:,:,2:15,ensembleMember); % 192 x 32 x 15
    
    %% Compare data to model
    
    % Field 1, SST Anomalies
    field = 1;
    diff_Annual_temp =  SST_Annual_Anomalies - model_SST_Annual_Anomalies(:,:,ensembleMember); % surface level
    diff(:,:,1,field,ensembleMember) = diff_Annual_temp ; % surface level, 192 x 32
    
    
    % Field 2, d18O calcite anomalie, planktonic
    field = 2;
    diff_d18O_calciteAnomalie_Planktonic = d18O_calciteAnomalie_Planktonic - model_d18O_calciteAnomalie(:,:,1,ensembleMember);
    diff(:,:,1,field,ensembleMember) = diff_d18O_calciteAnomalie_Planktonic; % 192 x 32
    
    
    % Field 3, d18O calcite anomalie, benthic
    field = 3;
    diff_d18O_calciteAnomalie_Benthic = d18O_calciteAnomalie_Benthic(:,:,2:15) - model_d18O_calciteAnomalie(:,:,2:15,ensembleMember);
    diff(:,:,2:15,field,ensembleMember) = diff_d18O_calciteAnomalie_Benthic; % 192 x 32 x 14
    
    
end

    %% Compute obs error cov matrix
    
    SigmaSquareTemp = [reshape(SIGMA_Annual.^2,[192*32,1]);  reshape(SIGMA_calciteAnomalie_Planktonic.^2,[192*32,1]);  reshape(SIGMA_calciteAnomalie_Benthic.^2,[192*32*15,1])];
    SigmaSquareTemp(isnan(SigmaSquareTemp)) = [];

    R = diag(SigmaSquareTemp);

end

