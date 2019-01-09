function Dimensions = sizeDimensions(config)


%% Check input arguments
narginchk(1,1);
if ischar(config.DimensionOrder) % field is a single string
    DimensionOrder = strsplit(config.DimensionOrder, ',');
else % field is a cell array of strings
    DimensionOrder = config.DimensionOrder;
end


%% Determine size of each dimension
nDimensions = length(DimensionOrder);
Dimensions = zeros(1,nDimensions);
for index = 1:nDimensions
    if isfield(config, DimensionOrder{index})
        Dimensions(index) = config.(DimensionOrder{index});
    else
        error('Dimension name %s not found as field in config file', DimensionOrder{index});
    end
end