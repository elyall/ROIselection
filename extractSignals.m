function [ROIdata, Data, Neuropil, ROIindex] = extractSignals(Images, ROIdata, ROIindex, varargin)
% [ROIdata, Data, Neuropil, ROIindex] = extractSignals(Images, ROIdata, ROIid, varargin)
% INPUTS:
% Images - images or filename(s) (cell array of strings or string)
% ROIdata - ROIdata (struct) or filename (string)
% ROIid - indices of ROIs within ROIdata to average or 'all' or 'new'
% ARGUMENTS:
% 'GPU' - performs computations on the GPU
% 'save' - saves output 'ROIdata' struct to file
% 'saveFile' - follow with the filename of the file to save to (defaults to
% ROIFile input if a filename is input for second input)
% 'loadType' - follow with 'MemMap' or 'Direct' to specify how to load the
% image files when a filename is input as the first input
% 'MotionCorrect' - follow with 'MCdata' structure, filename of mat file to
% load 'MCdata' structure from, or true to prompt for file selection
% 'Frames' - follow with vector specifying indices of frames to analyze.
% Last value can be 'inf', specifying to analyze all frames after the
% previous designated frame (ex: default is [1, inf] specifying all
% frames).

Mode = 'Cell';              % 'GPU', 'Cell', 'Sparse'
loadType = 'Direct';        % 'MemMap' or 'Direct'
saveOut = false;            % true or false
saveFile = '';              % filename to save ROIdata output to (defaults to ROIFile if one is input)
MotionCorrect = false;      % false, filename to load MCdata from, or true to prompt for file selection
Channel = 1;                % channel to extract data from
Depth = 1;                  % depth to extract data from
FrameIndex = [1, inf];      % vector of frame indices
borderLims = [0,0,34,34];   % number of pixels to remove from edges when computing ROI means (top, bottom, left, right)

% Memory settings
portionOfMemory = 0.08;     % find 10% or less works best
sizeRAM = 32000000000;      % amount of memory on your computer (UNIX-only)

verbose = false;
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

% Images
if iscellstr(Images) % filename input
    ImageFiles = Images;
    ROIdata.files = ImageFiles;
    switch loadType
        case 'MemMap'
            [Images, loadObj] = load2P(ImageFiles, 'Type', 'MemMap', 'Images', 'all');
            if strcmp(loadObj.files(1).ext, '.sbx')
                Images = intmax(loadObj.Precision) - Images;
            end
            numFramesPerLoad = loadObj.Frames;
        case 'Direct'
            [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2, 'Double');
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
    if loadObj.Depth == 1
        totalFrames = sum([loadObj.files(:).Frames]);
    else
        totalFrames = sum(floor([loadObj.files(:).Frames]/loadObj.Depth)+(rem([loadObj.files(:).Frames],loadObj.Depth)>=Depth));
    end
else % numeric array input
    loadType = false;
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
end


%% Cycle through each ROI averaging pixels within each frame

% Initialize output
Data = nan(numROIs, totalFrames);
Neuropil = nan(numROIs, totalFrames);

% Cycle through frames computing average fluorescence
fprintf('Extracting signals for %d ROI(s) from %d frame(s): %s\n', numROIs, numFrames, ROIFile)
if verbose
    parfor_progress(numFrames);
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
                [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', currentFrames, 'Channel', Channel, 'Depth', Depth);
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
            
            fprintf('Finished frame: ');
            for findex = 1:numImages
                
                % Calculate fluorescence signal
                Data(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, DataMasks, gpuArray(Images(:,findex))), 1));
                Neuropil(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, NeuropilMasks, gpuArray(Images(:,findex))), 1));
                
                if verbose
                    parfor_progress; % Update status
                end
                
            end %findex
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
            [img, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', findex, 'Channel', Channel, 'Depth', Depth, 'Verbose', false, 'double'); %direct
            
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
                parfor_progress; % Update status
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
            [img, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', findex, 'Channel', Channel, 'Depth', Depth, 'Verbose', false, 'double'); %direct
            
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
                parfor_progress; % Update status
            end
            
        end
        
end %Mode
if verbose
    parfor_progress(0);
end
fprintf('\nFinished extracting signals for %d ROI(s) from %d frame(s)\nSession took: %.1f minutes\n', numROIs, numFrames, toc/60)


%% Distribute data to structure
for rindex = 1:numROIs
    if ~isfield(ROIdata.rois(1),'rawdata');
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
