function latMeans_csFields = LatLegendreFit2csField(legendreFit,yC,xC,ylat,latMeans)
% Written by C.Breitkreuz (last modified 31.01.2019)
% Function LatLegendreFit2csField gets a legendreFit of a function (i.e.,
% the coefficients), then
% - spreads the legendreFit vector over a lat lon field
% - interpolates lat lon field to cs grid (yC,xC)
% - sets nan values at very high latitudes to closeest value in latMeans

% Input:
% legendreFit       - legendre fit coefficients of latMeans 
% yC, xC            - grid coordinates, center of grid cell
% ylat              - latitude coordinates
% latMeans          - original latudional mean values that are approximated by
%                     legendre Fit, this is just needed to put in values at
%                     nan-locations at very high latitude after
%                     interpolation from lat-lon grid

% Output:
% latMeans_csFields - latitudinal means spread on cubed sphere field

xi= -179:2:180; % longitude coordinates

nlong = length(xi);
nlat = length(legendreFit);

% lat-long grid 
yiField = repmat(ylat,[nlong,1]);
xiField = repmat(xi',[1,nlat]);


% Expand over field
legendreFit_Fields = repmat(legendreFit',[nlong,1] );


%%  Field of latitudinal means back to cubed sphere grid

latMeans_csFields = interp2( yiField, xiField, legendreFit_Fields, yC, xC);

% Get rid of nans;

% At his point NaNs are coming in, because the lat-lon-field doesn't
% have enough values at very high/low latitudes, i.e. lat = +-
% 87.940663871962499...

% I just fill these with the latitudional mean value at highest
% latitude, i.e. ylat(+-87)

temp = latMeans_csFields;
nanIndex = find(isnan(temp));
nanLat = yC(nanIndex);

% for Indian Icean, everything north of nanIndex gets set to value at
% nanIndex-1 (ylat(end))

for j = 1:length(nanLat)
    if nanLat(j) < ylat(1) % = -87
        temp(nanIndex(j)) = latMeans(1);
        % because ylat(1) = -87
    elseif nanLat(j) > ylat(end) % = 87
        temp(nanIndex(j)) = latMeans(nlat);
        % because ylat(59) = 87
    else
        error('Hier sollte kein Nan stehen!')
    end
end

latMeans_csFields(:,:) = temp;


end

