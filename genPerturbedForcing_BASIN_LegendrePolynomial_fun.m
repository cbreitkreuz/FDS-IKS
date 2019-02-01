function [EnsFields,parametersOut,latMeans,legendreFit,input_fields,legendreCoeffEns,ylat,Aqh_h2oField] ...
    = genPerturbedForcing_BASIN_LegendrePolynomial_fun(parametersIn, coeffIndeces, precipIndex, windIndex, basinIndicator,numParam, N, write2File,OutputPath,inputFileNames, inputFileNamesShort,smoothAlongBasins)

% genPerturbedForcing_BASIN_LegendrePolynomial_fun
% Written by C.Breitkreuz (last modified 31.01.2019)
% Generate ensemble of forcing fields for different atmospheric forcings
% from parameters (coefficients of a Legendre fit of the latitudional
% mean):
% - read original fields
% - compute latitudinal means (latMeans) (complicated for cubed sphere
% - compute Legendre fit of latudional means and interpolate back to cubed sphere


% Choose different (global or all three) basins for each parameter/field


% Input: (see Main_FDS_IKS.m)
% parametersIn - values for coefficient perturbations
% coeffIndeces - which coefficients for which field 
% precipIndex - indicates if and which field is precipitation, to delete negative values
% windIndex - indicates if and which field is wind fields, no legendre stuff implemented yet
% basinIndicator - array with 1s and 0s at size 1xnumFields, defines if parameters for field are to be estimated for individual basins (1) or globally (0)
% numParam - number of ctrl parameters per field and per basin! Total number of ctrl parameters is sum(numParam)
% N - order of Legendre polynomial/approximation - 1
% write2File - 1 or 0
% OutputPath - path for next iteration simulation 
% inputFileNames - original fields
% inputFileNamesShort - for naming new fields with perturbations
% smoothAlongBasins - smooth borders of ocean basins (e.g., between Atlantic and Indian Ocean) (1) or not (0)


% Output: (The outputs are not required in "Main_FDS_IKS.m" - just for
% convenience)
% EnsFields
% parametersOut - same as parametersIn just reshaped to maxNumParam x numberOfBasins x numFields x AllEnsSize
% latMeans
% legendreFit
% input_fields
% legendreCoeffEns
% ylat
% Aqh_h2oField


%% Reshape parameters

[numFields,numberOfBasins,maxNumParam] = size(coeffIndeces);

[~,AllEnsSize] = size(parametersIn);

parametersOut = nan(maxNumParam,numberOfBasins,numFields,AllEnsSize);

count = 0;

% reshape to maxNumParam x numberOfBasins x numFields x AllEnsSize
for i = 1:numFields
    for basin = 1:numberOfBasins
        for j = 1:numParam(i,basin)

            count = count + 1;
            parametersOut(j,basin,i,:) = parametersIn(count,:) ;
        end
    end
end


%% Read model grid

gridDirectory = <SPECIFY YOUR PATH HERE>

grids = rdmnc(fullfile(gridDirectory, 'grid.*'));

rA = grids.rA;
xC = grids.XC;
yC = grids.YC;

Nx = size(xC, 1);
Ny = size(xC, 2);

%% Read original forcing fields

inputDirectory = <SPECIFY YOUR PATH HERE>

% containing input fields for all different atmospheric forcings,i.e.
% numFields
input_fields = zeros(192,32,12,numFields);

for i = 1:numFields
    
    inputFile = fullfile(inputDirectory, inputFileNames{i});
    input_fields(:,:,:,i) = mit_readfield(inputFile, [Nx Ny 12], 'real*8');
end

size(input_fields)

%% Load extra fields if necessary (needs to be done for isotopic composition fields, see below)

UseRatio = zeros(numFields,1);

for i = 1:numFields
    
    
    % ratio precip
    if ~isempty( strfind(inputFileNames{i},'h2o18') ) % this really depends on name of field
        
        fprintf(['Compute precip ratio for field ',num2str(i),': ',inputFileNames{i},'\n'])
        UseRatio(i) = 1;
        inputFile = fullfile(inputDirectory, 'precip_h2o_LGM.bin');
        Precip_h2oField = mit_readfield(inputFile, [192 32 12], 'real*8');
        
        % compute ratio for precip
        input_fields(:,:,:,i)  = input_fields(:,:,:,i) ./ Precip_h2oField;
        
        % I need to take out very small values for precip h2o to make it work.
        % Otherwise the ratio h2o18/h2o has super high values and the balance
        % pkg blows everything up.
        errorIndex = find(Precip_h2oField < 10^-6);
        temp = input_fields(:,:,:,i);
        temp(errorIndex) = 2*10^-3;
        input_fields(:,:,:,i) = temp;
        
    end
    
    
    % ratio aqh
    if ~isempty( strfind(inputFileNames{i},'h218')  ) % this really depends on name of field
        
        fprintf(['Compute aqh ratio for field ',num2str(i),': ',inputFileNames{i},'\n'])
        UseRatio(i) = 1;
        
        inputFile = fullfile(inputDirectory, 'camiso_LGM_hum_h2o.bin');
        Aqh_h2oField = mit_readfield(inputFile, [192 32 12], 'real*8');
        
        aqhRatioConstant =  20/18 .* 2005.2 .* 10^-6;
        % compute ratio for aqh
        input_fields(:,:,:,i)  = ( input_fields(:,:,:,i) .* aqhRatioConstant )./ Aqh_h2oField;
        
    end
    
    % wind
    if ~isempty( strfind(inputFileNames{i},'wind')  )
        
        fprintf(['Also read vwind for field ',num2str(i),': ',inputFileNames{i},'\n'])
        
        inputFile = fullfile(inputDirectory, 'ccsm_vwind_LGM.bin');
        vwind = mit_readfield(inputFile, [192 32 12], 'real*8');
        vwind_basins = repmat(vwind,[1,1,1,numberOfBasins]);
        
    end

end


%% Read basin masks

% load mskBasU and mskBasV
MOC_prep = load('MOC_prep.mat');

% basin masks
mskBasC = MOC_prep.mskBasC;

mskBasC(mskBasC == 0) = nan;

if numberOfBasins == 3
    mskBasCAtlantic = mskBasC(:,:,1); % mask for Atlatic Ocean
    mskBasCIndic = mskBasC(:,:,2); % mask for Indian Ocean
    mskBasCPacific = mskBasC(:,:,3); % mask for Pacific Ocean
else
    mskBasCAtlantic = mskBasC(:,:,1); % mask for Atlatic Ocean
    mskBasCIndic = mskBasC(:,:,2); % mask for Indian Ocean
    mskBasCPacific = mskBasC(:,:,3); % mask for Pacific Ocean
    
    mskBasCPacificIndic = nan(192,32); % maks for Pacific + Indic
    mskBasCPacificIndic(mskBasCIndic==1) = 1;
    mskBasCPacificIndic(mskBasCPacific==1) = 1;
end

input_fields_basins = repmat(input_fields,[1,1,1,1,numberOfBasins]);


% add Mediterranean to Atlantic mask
index1 = find(xC < 39);
index2 = find(xC >-20);

index3 = find(yC < 45);
index4 = find(yC > 27);

index5 = intersect(index1,index2);
index6 = intersect(index3,index4);
MeditIndex = intersect(index5, index6);

mskBasCAtlantic(MeditIndex) = 1;


for i = 1:numFields
    %only the basin specific ones, for all other 1,2,3 are the same
    if basinIndicator(i)
        
        % split the fields up into three basins
        if numberOfBasins == 3
            input_fields_basins(:,:,:,i,1) = input_fields_basins(:,:,:,i,1) .* repmat(mskBasCAtlantic,[1 1 12]);
            input_fields_basins(:,:,:,i,2) = input_fields_basins(:,:,:,i,2) .* repmat(mskBasCIndic,[1 1 12]);
            input_fields_basins(:,:,:,i,3) = input_fields_basins(:,:,:,i,3) .* repmat(mskBasCPacific,[1 1 12]);
            
            if i == windIndex
                vwind_basins(:,:,:,1) = vwind_basins(:,:,:,1) .* repmat(mskBasCAtlantic,[1 1 12]);
                vwind_basins(:,:,:,2) = vwind_basins(:,:,:,2) .* repmat(mskBasCIndic,[1 1 12]);
                vwind_basins(:,:,:,3) = vwind_basins(:,:,:,3) .* repmat(mskBasCPacific,[1 1 12]);
            end
            
        elseif numberOfBasins == 2
            input_fields_basins(:,:,:,i,1) = input_fields_basins(:,:,:,i,1) .* repmat(mskBasCAtlantic,[1 1 12]);
            input_fields_basins(:,:,:,i,2) = input_fields_basins(:,:,:,i,2) .* repmat(mskBasCPacificIndic,[1 1 12]);
            
            if i == windIndex
                vwind_basins(:,:,:,1) = vwind_basins(:,:,:,1) .* repmat(mskBasCAtlantic,[1 1 12]);
                vwind_basins(:,:,:,2) = vwind_basins(:,:,:,2) .* repmat(mskBasCPacificIndic,[1 1 12]);
            end
            
        else
            error('Other number of basins not implemented.')
            
        end
    end
end


% save all masks in one variable for later (smoothing at edges of basins)
allMasks = zeros(192,32,numberOfBasins);

if numberOfBasins == 3
    allMasks(:,:,1) = mskBasCAtlantic;
    allMasks(:,:,2) = mskBasCIndic;
    allMasks(:,:,3) = mskBasCPacific;
elseif numberOfBasins == 2
    allMasks(:,:,1) = mskBasCAtlantic;
    allMasks(:,:,2) = mskBasCPacificIndic;
    
else
    error('Other number of basins not implemented.')
end

allMasks(isnan(allMasks)) = 0;


%% % % % %
%% Use broken lines to get yearly means along latitudinal bands (this is only for cubed sphere grid!!)
%% % % % %

%% Read broken lines file


% read broken lines file
bk_line_file = 'isoLat_cs32_59';
bkLine = load(bk_line_file);

nlat = length(bkLine.bkl_Ylat);
ylat = bkLine.bkl_Ylat;

% reshape into vector
nPg = Nx * Ny;
rA = reshape(rA,[nPg 1]);

input_fields_basins_yearlyMean = squeeze(nanmean(input_fields_basins,3)); % mean over months
input_vector = reshape(input_fields_basins_yearlyMean ,[nPg numFields numberOfBasins]);

% input_vector is size nPg x numFields x numberOfBasins, but for all fields where I dont
% use basins for the last dimension only the first is used. All three are
% the same for fields where I dont use basin.


if windIndex ~= 0
    vwind_basins_yearlyMean = squeeze(nanmean(vwind_basins,3)); % mean over months
    vwind_vector = reshape(vwind_basins_yearlyMean ,[nPg numberOfBasins]);
    
end

%% Build mean along broken lines for each month and each input field -> latitudinal means

% contains latitudinal means for each latitude
latMeans = zeros(nlat,numFields,numberOfBasins);
latMeans_vwind = zeros(nlat,numberOfBasins);

for lat = 1:nlat
    
    % 1,2,...,number of points for this latitude
    % Npts gives only e.g. 1:8, until -192 in bkl_IJuv
    Npts = 1:bkLine.bkl_Npts(lat);
    
    grid_indeces = bkLine.bkl_IJuv(Npts,lat); % indeces on cubed-sphere grid, between 1 and 192*32
    
    for i = 1:numFields
        
            
            if basinIndicator(i) == 1
                if numberOfBasins == 3
                    % three basins
                    basinsTemp = [1 2 3];
                elseif numberOfBasins==2
                    basinsTemp = [1 2];
                end
            else
                basinsTemp = 1;
                
            end
            
            for basin = basinsTemp % I need to do this only once for fields where I dont use basins
                
                % take only grid_indeces on current basin!!!! -> take out
                % nans (grid cells in other basin)
                
                %nanIndex flags grid cells from other basins!
                input_vector_temp = input_vector(grid_indeces,i,basin);
                nanIndexOtherBasins = find(isnan(input_vector_temp));
                
                if i == windIndex
                    vwind_vector_temp = vwind_vector(grid_indeces,basin);
                end
                
                % check if there are any values
                if length(nanIndexOtherBasins) == length(input_vector_temp)
                    warning(['No data in this mask for field: ',num2str(i),...
                        ' and for basin: ', num2str(basin), ' for latitude: ', num2str(ylat(lat))])
                end
                
                
                % set to one, value will be erased by multiplication with
                % rA
                input_vector_temp(nanIndexOtherBasins) = 1;
                
                if i == windIndex
                    vwind_vector_temp(nanIndexOtherBasins) = 1;
                end
                
                rAtemp = rA(grid_indeces);
                rAtemp(nanIndexOtherBasins) = 0;
                
                
                % build area-weighted mean for each latitude
                
                latMeans(lat,i,basin) = (rAtemp' * input_vector_temp)/sum(rAtemp);
                
                if i == windIndex
                    latMeans_vwind(lat,basin) = (rAtemp' * vwind_vector_temp)/sum(rAtemp);
                end
            end
            
    end % field
    
end % latitude


%% Compute Legendre fit of latudional means and interpolate back to cubed sphere

nanIndex = zeros(numFields,numberOfBasins);

anomalies = zeros(192,32,12,numFields,numberOfBasins);
legendreCoeff = zeros(N+1,numFields,numberOfBasins);
legendreFit = zeros(nlat,numFields,numberOfBasins);
latMeans_csFields_basins = zeros(Nx,Ny,numFields,numberOfBasins);


anomalies_vwind = zeros(192,32,12,numberOfBasins);
legendreCoeff_vwind= zeros(N+1,numberOfBasins);
legendreFit_vwind= zeros(nlat,numberOfBasins);
latMeans_csFields_basins_vwind = zeros(Nx,Ny,numberOfBasins);

for i = 1:numFields
    

        if basinIndicator(i) == 1
            if numberOfBasins == 3
                % three basins
                basinsTemp = [1 2 3];
            elseif numberOfBasins==2
                basinsTemp = [1 2];
            end
        else
            basinsTemp = 1;
            
        end
        
        for basin =  basinsTemp % I need to do this only once for fields where I dont use basins
            
            % Find index of nan values, where no data is available in this
            % basin, e.g. north of Indian Ocean
            
            nanIndexTemp = find(isnan(latMeans(:,i,basin)),1)-1;
            
            if isempty(nanIndexTemp )
                nanIndexTemp  = nlat;
            end
            
            nanIndex(i,basin) = nanIndexTemp;
            
            
            % Approximate latitudinal means with Legendre polynomials, get legendreFit and coefficients
            % The function legendrefit can be downloaded from MathWorks
            % File Exchange.
            [legendreCoeff(:,i,basin), legendreFit(1:nanIndex(i,basin),i,basin)] = legendrefit(latMeans(1:nanIndex(i,basin),i,basin), N, 'inv');
            
            if i == windIndex
                [legendreCoeff_vwind(:,basin), legendreFit_vwind(1:nanIndex(i,basin),basin)] = legendrefit(latMeans_vwind(1:nanIndex(i,basin),basin), N, 'inv');
            end
            
            
            %% This is just to check
            
            % Check if legendreCoeff fits to legendreFit (get function that
            % belongs to legendreCoefficients)
            f = LegendrePolynomialFromCoefficients(legendreCoeff(:,i,basin), nanIndex(i,basin));
            
            diff = f - legendreFit(1:nanIndex(i,basin),i,basin);
            
            if diff > 10^-12
                error('Something is wrong.')
            end
            
            %%  Expand lat means over field, interpolate to cubed sphere grid and compute anomalies
            
            % Function LatLegendreFit2csField expands legendreFit over a
            % lat-long grid, then interpolates to cs-grid.

            latMeans_csFields_basins(:,:,i,basin) = LatLegendreFit2csField(squeeze(legendreFit(1:nanIndex(i,basin),i,basin)),...
                yC,xC,ylat(1:nanIndex(i,basin)),...
                squeeze(latMeans(1:nanIndex(i,basin),i,basin)));
            
            % compute spatially and monthly varying anomalies
            for month = 1:12
                anomalies(:,:,month,i,basin) = squeeze(input_fields_basins(:,:,month,i,basin)) - squeeze(latMeans_csFields_basins(:,:,i,basin));
            end
            
            % for wind fields
            if i == windIndex
                latMeans_csFields_basins_vwind(:,:,basin) = LatLegendreFit2csField(squeeze(legendreFit_vwind(1:nanIndex(i,basin),basin)), ...
                    yC,xC,ylat(1:nanIndex(i,basin)),...
                    squeeze(latMeans_vwind(1:nanIndex(i,basin),basin)));
                
                
                for month = 1:12
                    anomalies_vwind(:,:,month,basin) = squeeze(vwind_basins(:,:,month,basin)) - squeeze(latMeans_csFields_basins_vwind(:,:,basin)); % 192 x 32 x 12 x numFields x numBasins=3
                end
                
            end
        end

end


%% Modify coefficients according to input parameters

% same parameter pEns(i) for all 12 month
legendreCoeffEns = zeros(N+1,numFields,numberOfBasins,AllEnsSize);
legendreFitEns = zeros(nlat,numFields,numberOfBasins,AllEnsSize);
latMeanFieldsEns = zeros(192,32,numFields,numberOfBasins,AllEnsSize);
EnsFieldsBasins = zeros(192,32,12,numFields,numberOfBasins,AllEnsSize);

legendreCoeffEns_vwind= zeros(N+1,numberOfBasins,AllEnsSize);
legendreFitEns_vwind= zeros(nlat,numberOfBasins,AllEnsSize);
latMeanFieldsEns_vwind= zeros(192,32,numberOfBasins,AllEnsSize);
EnsvWindBasin = zeros(192,32,12,numberOfBasins,AllEnsSize);

for i = 1:AllEnsSize
    
    % save original coefficient for those that don't get changed
    legendreCoeffEns(:,:,:,i) = legendreCoeff; % N+1 x numFields x ensSize x numBasins=3
    
    if windIndex ~= 0
        legendreCoeffEns_vwind(:,:,i) = legendreCoeff_vwind; % N+1  x numBasins=3 x ensSize
    end
    
    for j = 1:numFields
        
            
            if basinIndicator(j) % basin specific
                
                for basin =  1:numberOfBasins
                    
                    % modify all parameters for this field and this basin
                    for param = 1:numParam(j,basin) 
                        
                        % Add ctrl parameters:
                        % parameters are perturbatuions of coefficients!
                        legendreCoeffEns(coeffIndeces(j,basin,param),j,basin,i) = legendreCoeff(coeffIndeces(j,basin,param),j,basin) + parametersOut(param,basin,j,i);
                        
                        if j == windIndex 
                            legendreCoeffEns_vwind(coeffIndeces(j,basin,param),basin,i) = legendreCoeff_vwind(coeffIndeces(j,basin,param),basin) + parametersOut(param,basin,j,i);
                        end
                        
                    end
                    
                    
                    % For each set of coefficients compute Legendre
                    % polynomial
                    legendreFitEns(1:nanIndex(j,basin),j,basin,i) = LegendrePolynomialFromCoefficients(legendreCoeffEns(:,j,basin,i),nanIndex(j,basin));
                    
                    % For each polynomial compute cs lat mean Field
                    latMeanFieldsEns(:,:,j,basin,i) = LatLegendreFit2csField(legendreFitEns(1:nanIndex(j,basin),j,basin,i),...
                        yC,xC,ylat(1:nanIndex(j,basin)),...
                        latMeans(1:nanIndex(j,basin),j,basin));
                    
                    % same for wind field
                    if j == windIndex
                        
                        legendreFitEns_vwind(1:nanIndex(j,basin),basin,i) = LegendrePolynomialFromCoefficients(legendreCoeffEns_vwind(:,basin,i),nanIndex(j,basin)); %nanIndex here instead of nlat = 59, because it can be 54 for Indian ocean
                        
                        latMeanFieldsEns_vwind(:,:,basin,i) = LatLegendreFit2csField(legendreFitEns_vwind(1:nanIndex(j,basin),basin,i),...
                            yC,xC,ylat(1:nanIndex(j,basin)),...
                            latMeans_vwind(1:nanIndex(j,basin),basin));
                    end
                    
                end
                
            else % not basin specifc, I dont need nanIndex because globally values all latitudes
                

                basin = 1;
                
                for param = 1:numParam(j,basin)
                    
                    % Add ctrl parameters:
                    % parameters are perturbatuions of coefficients!
                    legendreCoeffEns(coeffIndeces(j,basin,param),j,basin,i) = legendreCoeff(coeffIndeces(j,basin,param),j,basin)+ parametersOut(param,basin,j,i);
                    
                    if j == windIndex 
                        legendreCoeffEns_vwind(coeffIndeces(j,basin,param),basin,i) = legendreCoeff_vwind(coeffIndeces(j,basin,param),basin) + parametersOut(param,basin,j,i);
                    end
                end
                
                
                
                % For each set of coefficients compute Legendre polynomial
                legendreFitEns(:,j,basin,i) = LegendrePolynomialFromCoefficients(legendreCoeffEns(:,j,basin,i),nlat);
                
                % For each polynomial compute cs lat mean Field
                latMeanFieldsEns(:,:,j,basin,i) = LatLegendreFit2csField(legendreFitEns(:,j,basin,i),yC,xC,ylat,latMeans(:,j,basin));
                
                % same for wind field
                if j == windIndex 
                    legendreFitEns_vwind(:,basin,i) = LegendrePolynomialFromCoefficients(legendreCoeffEns_vwind(:,basin,i),nlat); 
                    
                    latMeanFieldsEns_vwind(:,:,basin,i) = LatLegendreFit2csField(legendreFitEns_vwind(:,basin,i),yC,xC,ylat, latMeans_vwind(:,basin));
                end
            end
            
            
            % Add anomalies to each spread-out lat means 
            for month = 1:12
                EnsFieldsBasins(:,:,month,j,:,i) = squeeze(latMeanFieldsEns(:,:,j,:,i)) + squeeze(anomalies(:,:,month,j,:));
                
            end
            
    end
    
    
    for month = 1:12
        
        if  windIndex~= 0
            EnsvWindBasin(:,:,month,:,i) =  squeeze(latMeanFieldsEns_vwind(:,:,:,i)) + squeeze(anomalies_vwind(:,:,month,:));
        end

    end
    
end



%% Put the basins (perturbed with different parameters) back together

EnsFields = zeros(192,32,12,numFields,AllEnsSize);
vwindEns = zeros(192,32,12,AllEnsSize);

if numberOfBasins == 3
    
    basinIndexAtlantic = find(mskBasCAtlantic==1);
    basinIndexIndic = find(mskBasCIndic==1);
    basinIndexPacific = find(mskBasCPacific==1);
    
elseif numberOfBasins == 2
    
    basinIndexAtlantic = find(mskBasCAtlantic==1);
    basinIndexPacificIndic = find(mskBasCPacificIndic==1);
    
end

for month = 1:12
    for j = 1:numFields
        if basinIndicator(j)
            
            % pre-fill (Mediterranean not included in any mask-basin)
            tempGlobal = input_fields(:,:,month,j);
            
            if j == windIndex
                tempGlobalVwind = vwind(:,:,month);
            end
            
            for i = 1:AllEnsSize
                
                if numberOfBasins == 3
                    
                    tempAtlantic = EnsFieldsBasins(:,:,month,j,1,i); % 192 x 32
                    tempIndic = EnsFieldsBasins(:,:,month,j,2,i);
                    tempPacific = EnsFieldsBasins(:,:,month,j,3,i);
                    
                    
                    tempGlobal(basinIndexAtlantic) = tempAtlantic(basinIndexAtlantic);
                    tempGlobal(basinIndexIndic) = tempIndic(basinIndexIndic);
                    tempGlobal(basinIndexPacific) = tempPacific(basinIndexPacific);
                    
                    if j == windIndex
                        
                        tempAtlanticVwind =    EnsvWindBasin(:,:,month,1,i);
                        tempIndicVwind = EnsvWindBasin(:,:,month,2,i);
                        tempPacificVwind = EnsvWindBasin(:,:,month,3,i);
                        
                        tempGlobalVwind(basinIndexAtlantic) = tempAtlanticVwind(basinIndexAtlantic);
                        tempGlobalVwind(basinIndexIndic) = tempIndicVwind(basinIndexIndic);
                        tempGlobalVwind(basinIndexPacific) = tempPacificVwind(basinIndexPacific);
                        
                    end
                    
                    
                elseif numberOfBasins == 2
                    
                    tempAtlantic = EnsFieldsBasins(:,:,month,j,1,i); % 192 x 32
                    tempPacificIndic = EnsFieldsBasins(:,:,month,j,2,i);
                    
                    tempGlobal(basinIndexAtlantic) = tempAtlantic(basinIndexAtlantic);
                    tempGlobal(basinIndexPacificIndic) = tempPacificIndic(basinIndexPacificIndic);
                    
                    if j == windIndex
                        tempAtlanticVwind =    EnsvWindBasin(:,:,month,1,i);
                        tempPacificIndicVwind = EnsvWindBasin(:,:,month,2,i);
                        
                        tempGlobalVwind(basinIndexAtlantic) = tempAtlanticVwind(basinIndexAtlantic);
                        tempGlobalVwind(basinIndexPacificIndic) = tempPacificIndicVwind(basinIndexPacificIndic);
                        
                    end
                    
                end
                EnsFields(:,:,month,j,i) = tempGlobal;
                
                if j == windIndex
                    vwindEns(:,:,month,i)= tempGlobalVwind;
                end
                
            end
           
        else
            % If nothing basin specific, just take first one.
            basin = 1;
            EnsFields(:,:,month,j,:) = EnsFieldsBasins(:,:,month,j,basin,:);
            
            if j == windIndex
                vwindEns(:,:,month,:) = EnsvWindBasin(:,:,month,basin,:);
            end
        end
    end
end


%% Smoothing along ocean basin borders

if smoothAlongBasins
    
    [EnsFields,vwindEns] = smoothBasinBorders(numberOfBasins,allMasks,grids,basinIndicator,windIndex,numParam,EnsFields,vwindEns);
    
end

%% Transform ratio back


for i = 1:numFields
    
    if UseRatio(i) == 1
        if ~isempty( strfind(inputFileNames{i},'h2o18') )
            
            fprintf(['Transform ratio back for field ',num2str(i),': ',inputFileNames{i},'\n'])
            
            for particle = 1:AllEnsSize
                % EnsField is actually h2o18/h2o, but will be saved in h2o18.
                % In the model it will be divided by h2o, so ratio =  h2o18/h2o/h2o, so
                % I have to multiply by h2o, such that in the model ratio =
                % ((h2o18/h2o) * h2o) / h2o.
                EnsFields(:,:,:,i,particle) =  EnsFields(:,:,:,i,particle)  .* Precip_h2oField;
                
            end
            
        elseif ~isempty( strfind(inputFileNames{i},'h218') )
            
            fprintf(['Transform ratio back for field ',num2str(i),': ',inputFileNames{i},' with aqh ratio. ','\n'])
            
            for particle = 1:AllEnsSize
                
                % EnsField is actually (h2o18*ratio)/h2o, but will be saved in h2o18.
                % In the model it will be multiplied by ratio/h2o, so ratio =  (h2o18*ratio)/h2o * (ratio/h2o), so
                % I have to multiply by h2o/ratio, such that in the model ratio =
                % (( (h2o18*ratio)/h2o) * (h2o/ratio)) * (ratio/ h2o) .
                EnsFields(:,:,:,i,particle) =  EnsFields(:,:,:,i,particle) .* (Aqh_h2oField./ aqhRatioConstant);
            end
            
        end
    end
    
    
end


%% Take out negative values for precip (this can happen when anomalies are added)

% Set values to zero for all negativ values

if precipIndex ~=0
    temp = EnsFields(:,:,:,precipIndex,:);
    
    temp(temp<0) = 0;
    EnsFields(:,:,:,precipIndex,:) = temp;
end

%% Write everything to files


if write2File
    
    fprintf('Writing files!')
    

    for i = 1:AllEnsSize
        
        for j = 1:numFields
            
            fileID = fopen(fullfile(OutputPath,[inputFileNamesShort{j},num2str(i),'.bin']), 'w', 'b');
            fwrite(fileID, EnsFields(:,:,:,j,i), 'real*8');
            fclose(fileID);
            
            if j == windIndex
                fileID = fopen(fullfile(OutputPath,['ccsm_vwind_LGM_',num2str(i),'.bin']), 'w', 'b');
                fwrite(fileID, vwindEns(:,:,:,i) , 'real*8');
                fclose(fileID);
            end
            
        end
        
        
    end
    
    cd(OutputPath)
    save('parameters.mat','parametersOut')
    % back to current run directory
    cd <SPECIFY YOUR PATH HERE>
    
end

end

