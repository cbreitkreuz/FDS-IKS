function [EnsFields,vwindEns] = smoothBasinBorders(numberOfBasins,allMasks,grids,basinIndicator,windIndex,numParam,EnsFields,vwindEns)


xC = grids.XC;  % x-coordinate center
yC = grids.YC;

totalNumParams = sum(numParam(:));

[Nx,Ny,~,numFields,AllEnsSize]=size(EnsFields);
ensSize = (AllEnsSize -1)./totalNumParams;

% Use masks to find grid indeces that are at the border of a basin
% (basinborderIndex)

basinBorderIndex= zeros(400,numberOfBasins);

for basin = 1:numberOfBasins
    count = 1;
    
    % % %%% % % % % % % %%% % % % %% % %%% % % % %% % %%% % % % %% % %%% % % % %% % %%% % % % %
    % find all grid cells that have neighbour that is very
    % different from them in fieldTemp
    % % %%% % % % %% % %%% % % % %% % %%% % % % %% % %%% % % % %% % %%% % % % %
    
    maskTemp = allMasks(:,:,basin);
    
    %     plot_cs(maskTemp , grids, 'equidistant cylindrical','title');
    
    % find neighbour
    for i = 1:Nx*Ny
        
        % distances to all grid points
        dphi = r_dist(xC(i),yC(i),xC(:),yC(:));
        
        [sortdphi, sortIndex] = sort(dphi);
        
        % take 8 neighbours
        minIndex = sortIndex(2:9); % first one is grid cell itself
        %         minDist = sortdphi(2:9);
        
        % compare grid cell to 8 neighbours
        diffTemp = zeros(8,1);
        for k = 1:8
            
            diffTemp(k) = abs(maskTemp(i) - maskTemp(minIndex(k)));
            
        end
        
        % If there is a neighbour with big difference, save
        % index.
        if ~isempty(find(diffTemp > 0,1))
            basinBorderIndex(count,basin) = i;
            count = count +1;
            
            % add all of the neighbours as well
            basinBorderIndex(count:count+7,basin) = minIndex;
            count = count + 8;
        end
        
    end
    
end


%% for all basin border grid cells apply smoothing


% How many of the closest neighbours are used for smoothing
howManyNeighbours = 9;
% howManyNeighbours = 25;


for indexField = 1:numFields
    
    if basinIndicator(indexField)  % only do this for basin specific parameters
        
        
        for basin = 1:numberOfBasins
            
            
            thisBasinBorder = basinBorderIndex(:,basin);
            thisBasinBorder(thisBasinBorder==0) = [];
            
            %test = zeros(192,32);
            %test(thisBasinBorder) = 1;
            %plot_cs(test , grids, 'equidistant cylindrical','title');
            
            % find neighbours for this basin borders
            neighbourIndex = zeros(howManyNeighbours,length(thisBasinBorder));
            weights = zeros(howManyNeighbours,length(thisBasinBorder));
            
            for ii = 1:length(thisBasinBorder)
                
                % find neighbours for smoothing
                dphi = r_dist(xC(thisBasinBorder(ii)),yC(thisBasinBorder(ii)),xC(:),yC(:));
                [sortdphi, sortIndex] = sort(dphi);
                
                % take howManyNeighbours neighbours for
                % smoothing
                neighbourIndex(:,ii) = sortIndex(1:howManyNeighbours); % first one is grid cell itself
                minDist = sortdphi(1:howManyNeighbours);
                minDist(1) = 200; % write over distance 0 for next step. This is distance in kilometers, one gridcell ~ 280 km wide.
                
                % inverse distance weighting
                weights(:,ii) = 1./minDist;
                weights(:,ii) = weights(:,ii) ./sum(weights(:,ii)); % normalize
                
            end
            
            
            % find all particles that belong to this field and this
            % basin
            if AllEnsSize == 1 % this is the case if I don't put ensemble in but only one set of parameters
                particles = 1;
            else
                temp = numParam(1:indexField-1,:);
                StartParticle  = 2 + ensSize * sum(temp(:));  % from previous fields
                Add = (basin-1) * numParam(indexField,basin) * ensSize; % from previous basins
                particles = StartParticle+Add: StartParticle+Add + (ensSize * numParam(indexField,basin)) -1 ;
            end
            
            fprintf(['For field ',num2str(indexField),' and basin ',num2str(basin),', number of params: ',num2str(numParam(indexField,basin)), ' vary ensemble members ',num2str(particles),'\n'])
            
            
            for particle = particles
                for month = 1:12
                    
                    % take field that will be smoothed
                    fieldTemp = EnsFields(:,:,month,indexField,particle);
                    
                    if indexField == windIndex
                        vwind_fieldTemp = vwindEns(:,:,month,particle);
                    end
                    
                    % do smoothing for each grid cell at border
                    for ii = 1:length(thisBasinBorder)
                        
                        % do the smoothing
                        fieldTemp(thisBasinBorder(ii)) = weights(:,ii)' * fieldTemp(neighbourIndex(:,ii));
                        
                        if indexField == windIndex
                            vwind_fieldTemp(thisBasinBorder(ii)) = weights(:,ii)' * vwind_fieldTemp(neighbourIndex(:,ii));
                        end
                        
                        % put smoothed field back
                        EnsFields(:,:,month,indexField,particle) = fieldTemp;
                        
                        if indexField == windIndex
                            vwindEns(:,:,month,particle) = vwind_fieldTemp;
                        end
                        
                    end % border grid cells
                    
                end % month
            end % particles
            
        end % basins
    end
end



end

