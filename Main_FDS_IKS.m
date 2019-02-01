% Main_FDS_IKS
% Written by C.Breitkreuz 

% This skript
% - reads the model output of the current simulation
% - computes the model-data misfit
% - computes the new best guess parameters with the finite difference
%   sensitivity-iterative Kalman smoother (FDS-IKS)
% - computes the new forcing fields according to the state reduction
%   approach


% Below specify:
% - the paths to the model output and ctrl parameters from the last iteration
% - settings on ctrl parameters (example given)
% - paths to observational/proxy data
% - names of observational data files
% - your run diretory
% Additionally, specify grid directory, input directory and filenames in 
% genPerturbedForcing_BASIN_LegendrePolynomial_fun.m 
% and in getModelDataDiff_realLGM_SST_Anomalies_d18OAnomalies.m
% (All of these spots are highlighted with <SPECIFY YOUR PATH HERE>
% or <SPECIFY FILE NAME HERE>)



% Before running this skript you need :
% - to run the first model simulations (iteration 0)
% - the files: rdmnc (available in the MITgcm model code) and legendrefit
% (available via the Mathworks file exchange)

%% general settings

nextIteration = 10; % the last model simulation was "nextIteration - 1"
ExpName = 'LGM_Exp1';

write2File = 0; % 1 or 0

ensSize = 3; % number of sensitivty experiments per control parameter

smoothAlongBasins = 1; % smooth borders of ocean basins (e.g., between Atlantic and Indian Ocean) (1) or not (0)

%% settings about control parameters

% names of forcing files that control parameters belong to

inputFileNamesShort = {'ccsm_T_air_LGM_','ccsm_tot_precip_LGM_','precip_h2o18_LGM_','camiso_LGM_hum_h218o_'}; % for naming new fields with perturbations
inputFileNames  = {'ccsm_T_air_LGM.bin','ccsm_tot_precip_LGM.bin','precip_h2o18_LGM.bin','camiso_LGM_hum_h218o.bin'}; % original fields

numFields = length(inputFileNames); % number of different fields

precipIndex = 2; % indicates if and which field is precipitation, to delete negative values
windIndex = 0; % indicates if and which field is wind fields, no legendre stuff implemented yet

numberOfBasins = 2; % 3: Atlantic,Pacific and Indic OR 2: Atlantic, Pacific+Indic


basinIndicator = [1 1 1 1]; % array with 1s and 0s at size 1xnumFields, defines if parameters for field are to be estimated for individual basins (1) or globally (0)


% number of ctrl parameters per field and per basin! Total number of ctrl parameters is sum(numParam)
numParam = zeros(numFields,numberOfBasins);
numParam(1,1) = 3; % first field, first basin
numParam(1,2) = 3; % first field, second basin

numParam(2,1) = 1; % second field, first basin
numParam(2,2) = 1; %...

numParam(3,1) = 2;
numParam(3,2) = 2;

numParam(4,1) = 2;
numParam(4,2) = 2;


maxNumParam = max(numParam(:)); % max. number parameters for one field
totalNumParams = sum(numParam(:)); % total number of ctrl parameters

N = 2; % order of Legendre polynomial/approximation - 1

% coeffIndeces-1 is the index of the polynomial coefficient that we want to
% estimate, e.g. [1 3] -> coefficients c0 and c2, i.e. mean
% value and curvature of polynomial approximation.
% All other coefficients stay fixed.
coeffIndeces = nan(numFields,numberOfBasins,maxNumParam);

coeffIndeces(1,1,1) = 1; % first field
coeffIndeces(1,1,2) = 2;
coeffIndeces(1,1,3) = 3;

coeffIndeces(1,2,1) = 1;
coeffIndeces(1,2,2) = 2;
coeffIndeces(1,2,3) = 3;

coeffIndeces(2,1,1) = 1; % second field
coeffIndeces(2,2,1) = 1;

coeffIndeces(3,1,1) = 1;
coeffIndeces(3,1,2) = 3;

coeffIndeces(3,2,1) = 1; % third field
coeffIndeces(3,2,2) = 3;

coeffIndeces(4,1,1) = 1; % fourth field
coeffIndeces(4,1,2) = 3;

coeffIndeces(4,2,1) = 1;
coeffIndeces(4,2,2) = 3;


AllEnsSize = 1 + (totalNumParams * ensSize); % total number of simulations per iteration

% first guess of parameters (note that control parameters are actually the
% perturbation to the coefficients in the Legendre polynomial)
p0 = zeros(totalNumParams,1);


% First guess uncertainties of control parameters
% (this needs to be in opposite ordner with the dimensions such that the
% reshape below is doing the right thing)
sigma =nan(maxNumParam,numberOfBasins,numFields);

sigma(1,1,1)= 5;   % first field: air temp
sigma(2,1,1)= 5;
sigma(3,1,1)= 5;

sigma(1,2,1)= 5;
sigma(2,2,1)= 5;
sigma(3,2,1)= 5;

sigma(1,1,2)= 1*10^-9; % second field: precip
sigma(1,2,2)= 1*10^-9;

sigma(1,1,3)= 3*10^-6; % third field: isotopes precip
sigma(2,1,3)= 7*10^-6;

sigma(1,2,3)= 3*10^-6;
sigma(2,2,3)= 7*10^-6;

sigma(1,1,4)= 5*10^-6;  % fourth field: isotopes water vapor
sigma(2,1,4)= 9*10^-6;

sigma(1,2,4)= 5*10^-6;
sigma(2,2,4)= 9*10^-6;


sigma = reshape(sigma,[numFields*numberOfBasins*maxNumParam,1]);
sigma = sigma(~isnan(sigma));


% background error covariance matrix
P0 = zeros(totalNumParams,totalNumParams);

for i = 1:totalNumParams
    
    P0(i,i) = sigma(i).^2;
end


% Choose proxy data/observation files and files with their uncertainties:

%%%%%%%%%%%%%%%%%% Pseudo-proxy experiment

% targetFileOne = <SPECIFY FILE NAME HERE>
% targetFileOneWeights = 1/2^2; % 1/sigma^2, sigma = uncertainty of observational data 

% targetFileTwo = <SPECIFY FILE NAME HERE>
% targetFileTwoWeights = 1/0.22^2;  % 1/sigma^2, sigma = uncertainty of observational data 

% targetName = 'Target';

% path to target observations
% targetObsPath =  <SPECIFY YOUR PATH HERE>

%%%%%%%%%%%%%%%%%% LGM experiment

% do not need these for real LGM data
targetFileTwo = 'whatever';
targetFileTwoWeights = 'whatever';

targetFileOne= 'whatever';
targetFileOneWeights = 'whatever';

targetName = 'realLGM'; 

targetPathD18O = <SPECIFY YOUR PATH HERE>
targetPathSST = <SPECIFY YOUR PATH HERE>
%%%%%%%%%%%%%%%%%%

% Choose path to model output of last iteration and for next iteration
modelBasePath =[<SPECIFY YOUR PATH HERE>,ExpName,'/'];

OutputPath = [modelBasePath,ExpName,'_Iter',num2str(nextIteration)]; % for next iteration
modelPath = [modelBasePath,ExpName,'_Iter',num2str(nextIteration-1)]; % last iteration


%% Load parameters from previous iterations


% load parameters
fprintf(['Loading parameters from:', modelPath,'/parameters.mat','\n'])
load([modelPath,'/parameters.mat'])
parameters = parametersOut;


[maxNumParam2,numberBasins,numFields2,~]= size(parameters);
parameters = reshape(parameters,numFields2*maxNumParam2*numberBasins,AllEnsSize);

% take out nan rows
indexNan = find(isnan(parameters(:,1)));

parameters(indexNan,:) = [];
totalNumParam2 = size(parameters,1);


% just to check if parameters above correspond to the parameter.mat file
% that was read
if maxNumParam2 ~= maxNumParam || numFields2 ~= numFields || totalNumParam2 ~= totalNumParams
    error('MaxNumParam oder numFields wrong!!!!')
end


% perturbation experiment parameters from last iteration
param_Ens= reshape(parameters(:,2:AllEnsSize),[totalNumParams, ensSize, totalNumParams]);

% best guess parameters from last iteration
param_a = parameters(:,1);



%% Read model output

if isempty(strfind(targetName,'realLGM'))
    fprintf('Target is pseudo-proxy data! \n')
    
    [xFEnsTemp,DiffTemp,R]= getModelDataDiff(targetObsPath,targetFileOne,targetFileTwo,targetFileOneWeights,targetFileTwoWeights, modelPath,AllEnsSize);
    
    % xFEnsTemp - model output
    % DiffTemp - model-data difference
    % R - observation error covariance matrix
    
else
    fprintf('Target is real LGM data! \n')
    [xFEnsTemp,DiffTemp,R]= getModelDataDiff_realLGM_SST_Anomalies_d18OAnomalies(modelPath,AllEnsSize,targetPathSST,targetPathD18O);
end


% reshape to vector
[Nx,Ny,Nr,numObsFields,~] = size(xFEnsTemp);
xFEnsTemp2 = reshape(xFEnsTemp,[Nx*Ny*Nr*numObsFields,AllEnsSize]);
xF_DiffTemp = reshape(DiffTemp(:,:,:,:,1),[Nx*Ny*Nr*numObsFields,1]); % take only Diff from xF (simulation with best guess ctrl parameters)

% throw out nans
indexNoObs = find(isnan(xF_DiffTemp));
xF_DiffTemp(indexNoObs) = []; % take only those with observations
xF_Diff= xF_DiffTemp;

numberObs = length(xF_DiffTemp);
xFEnsTemp3 = zeros(numberObs,AllEnsSize);
for j = 1:AllEnsSize
    xFEnsTemp4 =  xFEnsTemp2(:,j);
    xFEnsTemp4(indexNoObs) = [];
    xFEnsTemp3(:,j) = xFEnsTemp4; % take only those with observations
end


% model output in observation space
xF = xFEnsTemp3(:,1); % model output simulation with best guess ctrl parameters
xFEns = reshape(xFEnsTemp3(:,2:end),[numberObs,ensSize,totalNumParams]); % model output perturbation experiments


%% Compute G, sensitivity matrix (Approximation of tangent linear model)

% AllEnsSize = (ensSize x totalnumParam)+1
% xFEns - numberObs x ensSize x numParam

G = zeros(numberObs,totalNumParams);

if ensSize == 1
    
    for i = 1:totalNumParams
        
        G(:,i) = (xFEns(:,1,i) - xF) ./ (param_Ens(i,:,i)- param_a(i));
    end
    
else % linear regression
    
    for i = 1:totalNumParams
        
        for j = 1:numberObs
            
            % Runs with perturbations
            % If ensSize = 3, use 4 simulations in regression ( + 1 base
            % line simulation)
            
            y(1) = xF(j);
            y(2:ensSize+1,1) = xFEns(j,:,i)';
            
            % perturbations
            x(1)  = param_a(i);
            x(2:ensSize+1,1)  = param_Ens(i,:,i)';
            
            X = [ones(length(x),1) x];
            
            b = X\y;
            
            G(j,i) = b(2); % slopes of linear regression
            
        end
    end
    
end


%% FDS-IKS

fprintf('Use FDS-IKS!')

K = P0 *  G' * inv(G * P0 * G' + R); % Kalman gain
param_a_New = p0 + K * (xF_Diff - G * (p0 - param_a) ); % new best guess parameters
Pa_New = (eye(totalNumParams) - K * G) * P0; % new error covariance matrix


%% Save stuff

if write2File
    
    % go to directory for next iteration
    cd(OutputPath)
    
    save('pA.mat','Pa_New')
    save('G.mat','G') % numberObs x totalNumberParams
    save('xF.mat','xF') %numberObs
    save('xF_Diff.mat','xF_Diff') % numberObs x 1
    save('xFEns.mat','xFEns') % numberObs x ensSize x totalNumParams
    
    % change back to this directory
    cd <SPECIFY YOUR PATH HERE>
    
end


%% Perturb new parameters

% start perturbation from new analysis

param_EnsNew = repmat(param_a_New,[1,ensSize,totalNumParams]);

perturb = zeros(ensSize,totalNumParams);
for i = 1:totalNumParams
    
    perturbScale = 100;
    
    perturb(:,i) = sigma(i)./perturbScale .* randn(ensSize,1); % compute perturbation
    
    param_EnsNew(i,:,i) = repmat(param_a_New(i),[ensSize, 1]) + perturb(:,i); % new parameters for sensitivity experiments
end


%% Generate new forcing fields from new parameters

newParameters = zeros(totalNumParams, AllEnsSize);

newParameters(:,1) = param_a_New;
newParameters(:,2:AllEnsSize) = reshape(param_EnsNew,[totalNumParams, ensSize * totalNumParams]);


genPerturbedForcing_BASIN_LegendrePolynomial_fun(newParameters, coeffIndeces, precipIndex, windIndex, basinIndicator,numParam, N, write2File, OutputPath, inputFileNames, inputFileNamesShort, smoothAlongBasins);

fprintf('Done!!! \n');

