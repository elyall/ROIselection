function [ROIMasks, NeuropilMasks, ROIdata] = createMasks(ROIMasks, varargin)

neuropiltype = 'dilation'; % 'sbx' or 'dilation'
ROIdata = [];

% Dilation only
NeuropilNum = 9;  % structuring element for creating neuropil masks' outer diameter during dilation step
BorderNum = 4;    % strucutring element for creating neuropil masks' inner diameter (border region)

% saving variables
ROIindex = [1 inf];
override = false;
saveOut = false;
saveFile = '';

verbose = false;
directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIdata'
                ROIdata = varargin{index+1};
                index = index + 2;
            case 'type'
                neuropiltype = varargin{index+1};
                index = index + 2;
            case 'NeuropilNum'
                NeuropilNum = varargin{index+1};
                index = index + 2;
            case 'BorderNum'
                BorderNum = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            case {'Verbose','verbose'}
                verbose = true;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('ROIMasks', 'var') || isempty(ROIMasks)
    [ROIMasks,p] = uigetfile({'*.rois;*.segment'},'Select ROI file:',directory);
    if isnumeric(ROIMasks)
        return
    end
    ROIMasks = fullfile(p, ROIMasks);
end


%% Load data
if ischar(ROIMasks)
    ROIFile = ROIMasks;
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
    [~,~,ext] = fileparts(ROIFile);
    switch ext
        case '.rois'
            load(ROIFile, 'ROIdata', '-mat');
            ROIMasks = reshape(full([ROIdata.rois(:).pixels]), size(ROIdata.rois(1).pixels,1), size(ROIdata.rois(1).pixels,2), numel(ROIdata.rois));
        case '.segment'
            load(ROIFile, 'mask', 'dim', '-mat');
            if issparse(mask)
                ROIMasks = reshape(full(mask), dim(1), dim(2), size(mask,2));
            else
                ROIMasks = mask;
            end
    end
elseif isstruct(ROIMasks)
    ROIdata = ROIMasks;
    ROIMasks = reshape(full([ROIdata.rois(:).pixels]), size(ROIdata.rois(1).pixels,1), size(ROIdata.rois(1).pixels,2), numel(ROIdata.rois));
end
[H, W, numROIs] = size(ROIMasks);


%% Determine ROIs to compute masks for
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(end-1)+1:size(ROIMasks,3));
end
if iscolumn(ROIindex)
    ROIindex = ROIindex';
end


%% Remove overlapping regions from ROI masks
countMatrix = sum(ROIMasks, 3);                    % determine number of ROIs found in each pixel
overlap = countMatrix > 1;                         % define regions where ROIs overlap
ROIMasks(repmat(overlap, 1, 1, numROIs)) = 0;      % remove regions of overlap from ROI masks


%% Create neuropil masks
fprintf('Creating ROI & Neuropil masks...');
countMatrix = logical(countMatrix);
NeuropilMasks = false(H,W,numROIs);
switch neuropiltype
    case 'dilation'
        
        numPixels = 2*prctile(sum(reshape(ROIMasks,H*W,numROIs)),95);
        if ~isempty(BorderNum)
            mask = imdilate(countMatrix, strel('disk',BorderNum));
        else
            mask = countmatrix;
        end
        if verbose; parfor_progress(numROIs); end;
        parfor rindex = 1:numROIs
            index = 0;
            while nnz(NeuropilMasks(:,:,rindex)) < numPixels
                NeuropilMasks(:,:,rindex) = imdilate(ROIMasks(:,:,rindex), strel('disk',NeuropilNum+index));
                NeuropilMasks(:,:,rindex) = NeuropilMasks(:,:,rindex).*~mask;
                index = index + 1;
            end
            if verbose; parfor_progress; end
        end
        if verbose; parfor_progress(0); end

    case 'sbx' % SBX method
        
        g = exp(-(-10:10).^2/2/2^2);
        maskb = conv2(g,g,double(countMatrix),'same')>.15; %.02                        % dilation for border region around ROIs (.15 -> 3 pixels)
        [xi,yi] = meshgrid(1:W,1:H);
        centroids = reshape([ROIdata.rois(ROIindex).centroid], 2, numROIs)';
        for rindex = 1:numROIs
            for neuropilrad = 40:5:100
                M = (xi-centroids(rindex,1)).^2+(yi-centroids(rindex,2)).^2 < neuropilrad^2;    % mask of pixels within the radius from the centroid (circle)
                NeuropilMasks(:,:,rindex) = M.*~maskb;                                          % remove ROIs and border regions
                if nnz(NeuropilMasks(:,:,rindex)) > 4000
                    break
                end
            end
        end

end %switch neuropiltype
fprintf('\tComplete\n');


%% Select a subsampling of neuropil pixels to ensure an equal signal to noise
% DiffinPixels = sum(reshape(NeuropilMasks, Height*Width, numROIs), 1) - sum(reshape(ROIMasks, Height*Width, numROIs), 1);
% for rindex = find(DiffinPixels > 0)
%     [I, J] = find(NeuropilMasks(:,:,rindex));
%     subset = datasample([I, J], DiffinPixels(rindex), 'Replace', false);
%     indices = sub2ind([Height, Width, numROIs], subset(:,1), subset(:,2), repmat(rindex, DiffinPixels(rindex), 1));
%     NeuropilMasks(indices) = 0;
% end
        

%% Distribute to structure
if nargout > 2 || (saveOut && ~isempty(saveFile))
    if isempty(ROIdata)
        load(saveFile, 'ROIdata', '-mat')
        if isempty(ROIdata)
            warning('Did not save output to file because initial ROIdata not given.');
            return
        end
    end
    
    % Distribute to structure
    for rindex = 1:numROIs
        if ~isfield(ROIdata.rois, 'mask') || isempty(ROIdata.rois(ROIindex(rindex)).mask) || ~isequal(ROIdata.rois(ROIindex(rindex)).mask, ROIMasks(:,:,rindex)) || override
            ROIdata.rois(ROIindex(rindex)).rawdata = [];
            ROIdata.rois(ROIindex(rindex)).mask = sparse(ROIMasks(:,:,rindex));
        end
        if ~isfield(ROIdata.rois, 'neuropilmask') || isempty(ROIdata.rois(ROIindex(rindex)).neuropilmask) || ~isequal(ROIdata.rois(ROIindex(rindex)).neuropilmask, NeuropilMasks(:,:,rindex)) || override
            ROIdata.rois(ROIindex(rindex)).rawneuropil = [];
            ROIdata.rois(ROIindex(rindex)).neuropilmask = sparse(NeuropilMasks(:,:,rindex));
        end
    end
end


%% Save output
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n',saveFile);
end

