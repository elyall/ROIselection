function [ROIdata, Data, Neuropil, ROIindex] = extractSignals(Images, ROIdata, ROIindex, varargin)
%EXTRACTSIGNALS Averages over all pixels within each ROI over time.
%   ROIdata = extractSignals() prompts user to select an images file, ROI
%   file, and averages over all ROIs for all time points. Returns ROIdata
%   with added sub-fields rois.rawdata and rois.rawneuropil.
%
%   ROIdata = extractSignals(IMAGES, ROIdata) averages all ROIs for all
%   time points in IMAGES. IMAGES can be a filename, a cell array of
%   strings specifying multiple files, or a matrix of dimension [H x W x D
%   x C x F]. ROIdata can be a string or an ROIdata struct.
%
%   ROIdata = extractSignals(...,ROIINDEX) averages over only ROIs
%   specified in ROIINDEX.
%
%   [ROIdata, DATA, NEUROPIL, ROIINDEX] = extractSignals(...) returns
%   matrices DATA and NEUROPIL of dimensions [numROIs x numFrames]
%   containing the rawdata and rawneuropil signals for all ROIs requested.
%   Also returns vector ROIINDEX of length numROIs which specifies the ROI
%   which specifies the ROI indices for each row of DATA and NEUROPIL.
%   
%   [...] = extractSignals(...,'verbose') displays status bar.
%
%   Other parameters can be specified via:
%   [...] = extractSignals(...,'PARAM1','VALUE1',...)
%   Parameters      Value 
%   'Mode'          'Cell', 'GPU', or 'SPARSE' specifying which method to 
%                   use to average over ROIs. (default = 'Cell')
%   'loadType'      'MemMap' or 'Direct' specifying how to load IMAGES.
%                   (default = 'Direct')
%   'MotionCorrect' MCdata struct, filename to load MCdata struct from,
%                   true to prompt for file selection, or false to not 
%                   perform motion correction. (default = false)
%   'Channel'       index of channel to average over. (default = 1)
%   'Depth'         index of depth to average over. (default = 1)
%   'Frames'        indices of frames to average for. (default = [1 inf])
%   'border'        vector of length 4 specifying number of pixels to
%                   ignore from [top, bottom, left, right]
%   'Save'          true to save ROIdata to file, false to not. (default =
%                   false)
%   'SaveFile'      filename of file to save ROIdata to. (default is 
%                   filename input if one is input)
%


Mode = 'Batch';              % 'GPU', 'Cell', 'Sparse'
loadType = 'Direct';        % 'MemMap' or 'Direct'
saveOut = false;            % true or false
saveFile = '';              % filename to save ROIdata output to (defaults to ROIFile if one is input)
MotionCorrect = false;      % false, filename to load MCdata from, or true to prompt for file selection
Channel = 1;                % channel to extract data from
Depth = 1;                  % depth to extract data from
FrameIndex = [1, inf];      % vector of relative frame indices for requested depth
borderLims = [0,0,4,0];     % number of pixels to ignore from edges when computing ROI means, inclusive (top, bottom, left, right)
verbose = false;            % booleon specifying whether to show status bar

% Memory settings
portionOfMemory = 0.08;     % find 10% or less works best
sizeRAM = 32000000000;      % amount of memory on your computer (UNIX-only)

% Placeholders
directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Mode'
                Mode = varargin{index+1};
                index = index + 2;
            case 'loadType'
                loadType = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            case 'MotionCorrect'
                MotionCorrect = varargin{index+1};
                index = index + 2;
            case 'Channel'
                Channel = varargin{index+1};
                index = index + 2;
            case 'Depth'
                Depth = varargin{index+1};
                index = index + 2;
            case {'Frames', 'frames', 'FrameIndex'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'borderLims','Border','border'}
                borderLims = varargin{index+1};
                index = index + 2;
            case {'Verbose', 'verbose'}
                if length(varargin)>index && islogical(varargin{index+1})
                    verbose = varargin{index+1};
                    index = index + 2;
                else
                    verbose = ~verbose;
                    index = index + 1;
                end
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Images', 'var') || isempty(Images)
    [Images,p] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select image file(s):',directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    elseif iscell(Images)
        Images = fullfile(p, Images);
    elseif ischar(Images)
        Images = {fullfile(p, Images)};
    end
elseif ischar(Images)
    Images = {Images};
end

if ~exist('ROIdata', 'var') || isempty(ROIdata)
    [ROIdata,p] = uigetfile({'*.rois;*.mat'},'Select ROI file:',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p, ROIdata);
end

if ~exist('ROIindex', 'var') || isempty(ROIindex)
    ROIindex = 'all'; % 'all' or 'new' or vector of indices
end

if isequal(MotionCorrect, true) % prompt for file selection
    [MotionCorrect, p] = uigetfile({'*.mat'},'Choose Experiment file to load MCdata from:',directory);
    if ischar(MotionCorrect)
        MotionCorrect = fullfile(p, MotionCorrect);
    end
end


%% Load in Data and determine dimensions

% ROIs
if ischar(ROIdata) % filename input
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
else
    ROIFile = 'file 1';
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end
if ~isfield(ROIdata,'Config')
    warning('ROIdata does not have image ''Config'' struct attached; assuming numDepths=1 & frameRate=15.49.');
    ROIdata.Config.Depth = 1;
    ROIdata.Config.FrameRate = 15.49;
    ROIdata.depth = 1;
end
if Depth~=ROIdata.depth
    warning('Changing Depth from %d to %d to match ROIdata.depth.',Depth,ROIdata.depth);
    Depth = ROIdata.depth;
end

% Images
if iscellstr(Images) % filename input
    ImageFiles = Images;
    ROIdata.files = ImageFiles;
    switch loadType
        case 'MemMap'
            [Images, loadObj,Config] = load2P(ImageFiles, 'Type', 'MemMap', 'Images', 'all');
            if strcmp(loadObj.files(1).ext, '.sbx')
                Images = intmax(loadObj.Precision) - Images;
            end
            numFramesPerLoad = loadObj.Frames;
        case 'Direct'
            [Images, loadObj,Config] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2, 'IndexType', 'relative', 'Depth', Depth, 'Double');
            sizeFrame = whos('Images');
            sizeFrame = sizeFrame.bytes;
            if ispc
                mem = memory;
                numFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
            else
                numFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
            end
    end
    Height = loadObj.Height;
    Width = loadObj.Width;
    if Config(1).Depth == 1 % assumes all files have the same depth
        totalFrames = sum([loadObj.files(:).Frames]);
    else
        totalFrames = sum(floor([loadObj.files(:).Frames]/Config(1).Depth)+(rem([loadObj.files(:).Frames],Config(1).Depth)>=Depth));
    end
else % numeric array input
    loadType = false;
    ImageFiles = {};
    [Height, Width, ~, ~, totalFrames] = size(Images);
    numFramesPerLoad = totalFrames;
end

% Create border
if any(borderLims)
    Border = ones(Height,Width);
    Border(1:borderLims(1),:) = nan;            % top
    Border(end-borderLims(2)+1:end,:) = nan;    % bottom
    Border(:,1:borderLims(3)) = nan;            % left
    Border(:,end-borderLims(4)+1:end) = nan;    % right
    borderLims = true;
else
    borderLims = false;
end


%% Determine ROIs to extract signals for
if ischar(ROIindex)
    switch ROIindex
        case {'all', 'All'}
            ROIindex = 1:numel(ROIdata.rois);
        case {'new', 'New'}
            ROIindex = find(arrayfun(@(x) (isempty(x.rawdata)), ROIdata.rois));
    end
elseif isnumeric(ROIindex) && ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

if numROIs == 0
    fprintf('No new ROIs to extract signals from: %s\n', ROIFile);
    Data = [];
    Neuropil = [];
    return
end

if any(arrayfun(@(x) isempty(x.mask),ROIdata.rois(ROIindex)))
    warning('ROI masks did not exist for %d ROIs, running ''createMasks''...',nnz(arrayfun(@(x) isempty(x.mask),ROIdata.rois(ROIindex))));
    [~,~,ROIdata] = createMasks(ROIdata);
end


%% Determine frames to process
if FrameIndex(end)==inf
    FrameIndex = cat(2, FrameIndex(1:end-1), FrameIndex(end-1)+1:totalFrames);
end
numFrames = numel(FrameIndex);


%% Load in motion correction information
if ischar(MotionCorrect) % load in MCdata structure
    load(MotionCorrect, 'MCdata', '-mat')
    if ~exist('MCdata', 'var')
        MotionCorrect = false;
    else
        MotionCorrect = true;
    end
elseif isstruct(MotionCorrect) % MCdata structure input
    MCdata = MotionCorrect;
    MotionCorrect = true;
else
    MCdata = [];
end


%% Cycle through each ROI averaging pixels within each frame

% Initialize output
Data = nan(numROIs, totalFrames);
Neuropil = nan(numROIs, totalFrames);

% Cycle through frames computing average fluorescence
fprintf('Extracting signals for %d ROI(s) from %d frame(s): %s\n', numROIs, numFrames, ROIFile)
if verbose
    h = parfor_progress(numFrames);
end
tic;
switch Mode
    case 'GPU'
        
        % Define masks
        DataMasks = gpuArray(double(reshape([ROIdata.rois(ROIindex).mask], [Height*Width, numROIs])));
        NeuropilMasks = gpuArray(double(reshape([ROIdata.rois(ROIindex).neuropilmask], [Height*Width, numROIs])));
        DataMasks(~DataMasks) = NaN; % turn logical 0s to NaNs for NANMEAN
        NeuropilMasks(~NeuropilMasks) = NaN; % turn logical 0s to NaNs for NANMEAN
        
        % Cycle through frames in batches
        for bindex = 1:numFramesPerLoad:numFrames % direct loading only -> load frames in batches
            lastframe = min(bindex+numFramesPerLoad-1, numFrames);
            currentFrames = FrameIndex(bindex:lastframe);
            
            % direct loading only -> load current batch
            if strcmp(loadType, 'Direct')
                if bindex ~= 1
                    fprintf('\n');
                end
                [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', currentFrames, 'IndexType', 'relative', 'Channel', Channel, 'Depth', Depth);
            end
            
            % Remove border pixels
            if borderLims
                Images = bsxfun(@times, Images, Border);
            end
            
            % Correct for motion
            if MotionCorrect
                fprintf('\b\tCorrecting motion...');
                Images = applyMotionCorrection(Images, MCdata, loadObj);
                fprintf('\tComplete\n');
            end
            
            % Reshape images
            numImages = size(Images, 5);
            Images = double(reshape(Images(:,:,1,1,:), Height*Width, numImages));
            
            % Calculate fluorescence signal
            for findex = 1:numImages
                Data(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, DataMasks, gpuArray(Images(:,findex))), 1));
                Neuropil(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, NeuropilMasks, gpuArray(Images(:,findex))), 1));
                if verbose
                    parfor_progress(h); % Update status
                end
            end %findex
        end %bindex
        
    case 'Batch'
        
        % Define masks
        DataMasks = cell(numROIs, 1);
        NeuropilMasks = cell(numROIs, 1);
        for rindex = 1:numROIs
            DataMasks{rindex} = find(ROIdata.rois(ROIindex(rindex)).mask);
            NeuropilMasks{rindex} = find(ROIdata.rois(ROIindex(rindex)).neuropilmask);
        end
        
        % Cycle through frames in batches
        for bindex = 1:numFramesPerLoad:numFrames % direct loading only -> load frames in batches
            lastframe = min(bindex+numFramesPerLoad-1, numFrames);
            currentFrames = FrameIndex(bindex:lastframe);
            
            % direct loading only -> load current batch
            if strcmp(loadType, 'Direct')
                if bindex ~= 1
                    fprintf('\n');
                end
                [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', currentFrames, 'IndexType', 'relative', 'Channel', Channel, 'Depth', Depth, 'double');
            end
            
            % Remove border pixels
            if borderLims
                Images = bsxfun(@times, Images, Border);
            end
            
            % Correct for motion
            if MotionCorrect
                fprintf('\b\tCorrecting motion...');
                Images = applyMotionCorrection(Images, MCdata, loadObj);
                fprintf('\tComplete\n');
            end
                        
            % Calculate fluorescence signal
            numImages = size(Images, 5);
            parfor rindex = 1:numROIs
                inds = DataMasks{rindex} + Height*Width*(0:numImages-1); % indices of ROI in all frames
                Data(rindex,currentFrames) = nanmean(reshape(Images(inds(:)),[numel(DataMasks{rindex}),numImages]));
                inds = NeuropilMasks{rindex} + Height*Width*(0:numImages-1);
                Neuropil(rindex,currentFrames) = nanmean(reshape(Images(inds(:)),[numel(NeuropilMasks{rindex}),numImages]));
            end %rindex

%             if verbose
%                 parfor f = 1:numImages
%                     parfor_progress(h); % update status
%                 end
%             end
            
        end %bindex
        
    case 'Cell'
        
        % Define masks
        DataMasks = cell(numROIs, 1);
        NeuropilMasks = cell(numROIs, 1);
        for rindex = 1:numROIs
            DataMasks{rindex} = find(ROIdata.rois(ROIindex(rindex)).mask);
            NeuropilMasks{rindex} = find(ROIdata.rois(ROIindex(rindex)).neuropilmask);
        end
        
        % Cycle through frames
        parfor findex = FrameIndex
            
            % Load Frame
%             if loadType % fails to properly avoid else in parfor loop
                [img, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', findex, 'IndexType', 'relative', 'Channel', Channel, 'Depth', Depth, 'double'); %direct
%             else
%                 loadObj = []; % FIX LATER
%                 img = Images(:,:,Depth,Channel,findex);
%             end
            
            % Remove border pixels
            if borderLims
                img = img.*Border;
            end
            
            % Correct for motion
            if MotionCorrect
                img = applyMotionCorrection(img, MCdata, loadObj);
            end            
            
            % Calculate fluorescence signal
            for rindex = 1:numROIs
                Data(rindex,findex) = nanmean(img(DataMasks{rindex}));
                Neuropil(rindex,findex) = nanmean(img(NeuropilMasks{rindex}));
                % Neuropil(rindex,findex) = trimmean(img(NeuropilMasks{rindex}), 10); % sbx method
            end
            
            if verbose
                parfor_progress(h); % Update status
            end
            
        end %findex
        
    case 'Sparse'
        
        % Initialize output
        DataMasks = reshape([ROIdata.rois(ROIindex).mask], [Height*Width, numROIs]);
        NeuropilMasks = reshape([ROIdata.rois(ROIindex).neuropilmask], [Height*Width, numROIs]);
        if ~issparse(DataMasks)
            DataMasks = sparse(DataMasks);
        end
        if ~issparse(NeuropilMasks)
            NeuropilMasks = sparse(NeuropilMasks);
        end
        numPixelsData = full(sum(DataMasks));
        numPixelsNeuropil = full(sum(NeuropilMasks));
        
        % Cycle through frames
        parfor findex = FrameIndex
            
            % Load Frame
            [img, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', findex, 'IndexType', 'relative', 'Channel', Channel, 'Depth', Depth, 'double'); %direct
            
            % Remove border pixels
            if borderLims
                img = img.*Border;
            end
            
            % Correct for motion
            if MotionCorrect
                img = applyMotionCorrection(img, MCdata, loadObj);
            end
            
            % Calculate fluorescence signal
            Data(:,findex) = full(nansum(bsxfun(@times, DataMasks, img)))./numPixelsData;
            Neuropil(:,findex) = full(nansum(bsxfun(@times, DataMasks, img)))./numPixelsNeuropil;
            
            if verbose
                parfor_progress(h); % Update status
            end
            
        end
        
end %Mode
if verbose
    parfor_progress(h,0);
end
fprintf('\tComplete.\tSession took: %.1f minutes\n', toc/60)


%% Distribute data to structure
for rindex = 1:numROIs
    if ~isfield(ROIdata.rois(1),'rawdata')
        ROIdata.rois(1).rawdata = [];
    end
    if isempty(ROIdata.rois(ROIindex(rindex)).rawdata) % replace whole vector
        ROIdata.rois(ROIindex(rindex)).rawdata = Data(rindex, :);
        ROIdata.rois(ROIindex(rindex)).rawneuropil = Neuropil(rindex, :);
    else % replace frames that were computed
        ROIdata.rois(ROIindex(rindex)).rawdata(FrameIndex) = Data(rindex, FrameIndex);
        ROIdata.rois(ROIindex(rindex)).rawneuropil(FrameIndex) = Neuropil(rindex, FrameIndex);
    end
end


%% Save data to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    if exist('ImageFiles', 'var')
        save(saveFile, 'ImageFiles', '-mat', '-append');
    end
    if isequal(ROIindex, 1:numel(ROIdata.rois))
        save(saveFile, 'Data', 'Neuropil', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end
