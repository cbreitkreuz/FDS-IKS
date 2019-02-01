function [model_output, diff, R] = getModelDataDiff(targetPath,targetFileOne,targetFileTwo,targetFileOneWeights,targetFileTwoWeights,modelPath,ensSize)
% Written by C.Breitkreuz (last modified 31.01.2019)
% getModelDataDiff reads model output and computes model-data misift

% Input:
% targetPath - path to observational data
% targetFileOne
% targetFileTwo
% targetFileOneWeights
% targetFileTwoWeights
% modelPath - path to model output
% ensSize - total number of ensemble members


% Output:
% model_output 
% diff - model-data misfif
% R - observation error covariance matrix

%% Read model grid

gridDirectory = '/m/wrk3/cbreitkreuz/MITgcm_exp/global_ocean.cs32x15_LGM/run_LGM_wiso_2/mnc_output_all';

grids = rdmnc(fullfile(gridDirectory, 'grid.*'));

hFacC = grids.HFacC;

% Determine array sizes
Nx = size(hFacC, 1);
Ny = size(hFacC, 2);
Nr = size(hFacC, 3);

mask = hFacC;
mask(mask > 0) =  1;
mask(mask==0) = nan;


%% If only observational data from Pacific+Indian Ocean or Atlantic Ocean is supposed to be used:

% % load mskBasU and mskBasV
% MOC_prep = load('MOC_prep.mat');
%
% % basin masks
% mskBasC = MOC_prep.mskBasC;
% maskAtlantic=repmat(mskBasC(:,:,1),[1 1 15]); % mask for Atlatic Ocean
%
%
% % % %% % %%%%%%%%%%%%%
% % PACIFIC+INDIC MASK
%
% maskNotAtlantic = maskAtlantic;
% maskNotAtlantic(maskNotAtlantic > 0) = 2;
%
% maskNotAtlantic(maskNotAtlantic == 0) = 1;
% maskNotAtlantic(maskNotAtlantic == 2) = nan;
%
% mask = mask.* maskNotAtlantic;
%
% % % %% % %%%%%%%%%%%%%
%
% % ATLANTIC MASK
% mask = mask.* maskAtlantic;
%
% % % %% % %%%%%%%%%%%%%

%% Read target data and weights

TargetField_1 = mit_readfield(fullfile(targetPath, targetFileOne), [192 32 15], 'real*8') .* mask;
TargetField_2 = mit_readfield(fullfile(targetPath, targetFileTwo), [192 32 15], 'real*8') .* mask;

% read weights
if ~(ischar(targetFileOneWeights) && ischar(targetFileTwoWeights))
    
    % scalar
    WeightsTarget_One = targetFileOneWeights ;
    WeightsTarget_Two = targetFileTwoWeights ;
    
else % char -> file
    
    fprintf('Reading weight file!')
    WeightsTarget_One= mit_readfield(fullfile(targetPath, targetFileOneWeights), [192 32 15], 'real*8') .* mask;
    WeightsTarget_Two = mit_readfield(fullfile(targetPath, targetFileTwoWeights), [192 32 15], 'real*8') .* mask;
    
end

TargetField_1(WeightsTarget_One == 0) = nan;
WeightsTarget_One(WeightsTarget_One==0)=nan;
TargetField_2(WeightsTarget_Two == 0) = nan;
WeightsTarget_Two(WeightsTarget_Two==0)=nan;

length(find(~isnan(TargetField_1)))
length(find(~isnan(WeightsTarget_One)))
length(find(~isnan(TargetField_2)))
length(find(~isnan(WeightsTarget_Two)))

%% Start loop over ensemble members

model_output = zeros(Nx,Ny,Nr,2,ensSize); %192x32x15 x 12 month x (temp & salt/d18O) x ensSize
diff = zeros(Nx,Ny,Nr,2,ensSize);


% Read vector of iterations and find last time slice

a  = rdmnc(fullfile(fullfile(modelPath, '/mnc_output_all_1'), 'dynDiag.*'), 'iter');
it = a.iter(end);


for ensembleMember = 1:ensSize
    
    fprintf(['Reading ensemble member ',num2str(ensembleMember),'\n'])
    
    %% Set directory name according to ensemble member
    
    modelExperiment = ['/mnc_output_all_',num2str(ensembleMember)];
    
    modelDirectory  = fullfile(modelPath, modelExperiment);
    
    %% Read model temperature
    
    if ~isempty(strfind(targetFileOne,'Theta'))
        
        % Read data structure from file
        inputField = rdmnc(fullfile(modelDirectory, 'dynDiag.*'),'THETA', it);
        THETA = inputField.THETA;
        
        model_output(:,:,:,1,ensembleMember) = THETA;
        
    else
        error('First one should be theta file!')
        
    end
    
    %% Read model d18O or salinity
    
    if ~isempty(strfind(targetFileTwo,'Salt'))
        
        % Read data structure from file
        inputField = rdmnc(fullfile(modelDirectory, 'dynDiag.*'),'SALT', it);
        SALT =  inputField.SALT;
        
        
        model_output(:,:,:,2,ensembleMember) = SALT;
        
    elseif ~isempty(strfind(targetFileTwo,'d18O'))
        
        
        % Read data structure from file
        inputField_ptracer = rdmnc(fullfile(modelDirectory, 'tracerDiag.*'),'TRAC01','TRAC02',it);
        
        H2O16 = inputField_ptracer.TRAC01;
        H2O18 = inputField_ptracer.TRAC02;
        
        
        % Calculate d180
        model_output(:,:,:,2,ensembleMember) = (((H2O18./H2O16)./(1/498.7))-1)*1000;
        
    else
        error('Neither salt nor d18O File...!')
        
    end
    
    
    %% Calculate difference and mask!
    
    %192x32x15 x (temp & salt) x ensSize
    
    diff(:,:,:,1,ensembleMember) = ( TargetField_1 - model_output(:,:,:,1,ensembleMember)) .* mask; % 192 x 32 x 15
    diff(:,:,:,2,ensembleMember) = ( TargetField_2 - model_output(:,:,:,2,ensembleMember)) .* mask;
    
    
end % ensemble


%% Compute observational error cov matrix

numberObs(1) = length(find(~isnan(diff(:,:,:,1,ensembleMember))));
numberObs(2) = length(find(~isnan(diff(:,:,:,2,ensembleMember))));


if ~(ischar(targetFileOneWeights) && ischar(targetFileTwoWeights)) % spatially constant uncertainty given
    
    SigmaSquareTemp = [ repmat(1./WeightsTarget_One,[numberObs(1),1]); repmat(1./WeightsTarget_Two,[numberObs(2),1]) ];
    R = diag(SigmaSquareTemp);
    
else % field of uncertainties given
    
    SigmaSquareTemp = [ reshape(1./WeightsTarget_One,[192*32*15,1]) ;  reshape(1./WeightsTarget_Two,[192*32*15,1]) ];
    SigmaSquareTemp(isnan(SigmaSquareTemp)) = [];
    
    R = diag(SigmaSquareTemp);
    
end


end

