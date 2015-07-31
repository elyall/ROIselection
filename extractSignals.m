function [ROIdata, Data, Neuropil, ROIid] = extractSignals(Images, ROIdata, ROIid, varargin)
% [ROIdata, Data, Neuropil] = extractSignals2(Images, ROIdata, ROIid, varargin)
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

GPU = true; % true or false
loadType = 'Direct'; % 'MemMap' or 'Direct'
saveOut = true; % true or false
saveFile = ''; % filename to save ROIdata output to (defaults to ROIFile if one is input)
MotionCorrect = false; % false, filename to load MCdata from, or true to prompt for file selection
FrameIndex = [1, inf]; % vector of frame indices

% Memory settings
portionOfMemory = 0.08; % find 10% or less works best
sizeRAM = 32000000000; % amount of memory on your computer (UNIX-only)

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            case 'GPU'
                GPU = true;
                index = index + 1;
            case 'loadType'
                loadType = varargin{index+1};
                index = index + 2;
            case 'MotionCorrect'
                MotionCorrect = varargin{index+1};
                index = index + 2;
            case 'Frames'
                FrameIndex = varargin{index+1};
                index = index + 2;
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
    directory = CanalSettings('DataDirectory');
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
    directory = CanalSettings('DataDirectory');
    [ROIdata,p] = uigetfile({'*.rois;*.mat'},'Select ROI file:',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p, ROIdata);
end

if ~exist('ROIid', 'var') || isempty(ROIid)
    ROIid = 'all'; % 'all' or 'new' or vector of indices
end

if isequal(MotionCorrect, true) % prompt for file selection
    directory = CanalSettings('ExperimentDirectory');
    [MotionCorrect, p] = uigetfile({'*.mat'},'Choose Experiment file to load MCdata from:',directory);
    if ischar(MotionCorrect)
        MotionCorrect = fullfile(p, MotionCorrect);
    end
end

tic
%% Load in Data and determine dimensions

% ROIs
if ischar(ROIdata) % filename input
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
else
    ROIFile = '';
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
    Depth = loadObj.Depth;
    numFrames = sum([loadObj.files(:).Frames]);
else % numeric array input
    loadType = false;
    [Height, Width, Depth, ~, numFrames] = size(Images);
    numFramesPerLoad = numFrames;
end
if Depth > 1
    error('Currently does not work for datasets with multiple z-planes');
end


%% Determine ROIs to extract signals for
if ischar(ROIid)
    switch ROIid
        case 'all'
            numROIs = numel(ROIdata.rois);
            ROIid = 1:numROIs;
        case 'new'
            ROIid = find(arrayfun(@(x) (isempty(x.rawdata)), ROIdata.rois));
            numROIs = numel(ROIid);
    end
end


%% Determine frames to process
if FrameIndex(end)==inf
    FrameIndex = cat(2, FrameIndex(1:end-1), FrameIndex(end-1)+1:numFrames);
end
totalFrames = numel(FrameIndex);


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

% Create masks (1s of pixels to average, NaNs of pixels to ignore)
if GPU
    DataMasks = gpuArray(double(reshape([ROIdata.rois(ROIid).mask], [Height*Width, numROIs])));
    NeuropilMasks = gpuArray(double(reshape([ROIdata.rois(ROIid).neuropilmask], [Height*Width, numROIs])));
else
    DataMasks = double(reshape([ROIdata.rois(ROIid).mask], [Height*Width, numROIs]));
    NeuropilMasks = double(reshape([ROIdata.rois(ROIid).neuropilmask], [Height*Width, numROIs]));
end
DataMasks(~DataMasks) = NaN; % turn logical 0s to NaNs for NANMEAN
NeuropilMasks(~NeuropilMasks) = NaN; % turn logical 0s to NaNs for NANMEAN

% Initialize output
Data = nan(numROIs, numFrames);
Neuropil = nan(numROIs, numFrames);

% Cycle through frames computing average fluorescence
fprintf('Extracting signals from %d frames for %d ROIs: %s\n', totalFrames, numROIs, ROIFile)
for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
    lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
    currentFrames = FrameIndex(bindex:lastframe);
    
    % direct loading only -> load current batch
    if strcmp(loadType, 'Direct')
        [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', currentFrames); %direct
    end
    
    % Correct for motion
    if MotionCorrect
        fprintf('Correcting motion...');
        Images = applyMotionCorrection(Images, MCdata, loadObj);
        fprintf('\tComplete\n');
    end
    
    fprintf('\tFinished frame:');
    for findex = 1:size(Images, 5);
                
        % Calculate fluorescence signal
        img = double(reshape(Images(:,:,1,1,findex), Height*Width, 1));
        if GPU
            Data(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, DataMasks, img), 1));
            Neuropil(:,currentFrames(findex)) = gather(nanmean(bsxfun(@times, NeuropilMasks, img), 1));
        else
            Data(:,currentFrames(findex)) = nanmean(bsxfun(@times, DataMasks, img), 1);
            Neuropil(:,currentFrames(findex)) = nanmean(bsxfun(@times, NeuropilMasks, img), 1);
        end
        
        % Update status
        if ~mod(currentFrames(findex), 20)
            fprintf('\t%d', currentFrames(findex))
            if ~mod(currentFrames(findex), 600)
                fprintf('\n\tFinished frame:');
            end
        end
        
    end %findex
end %bindex

for rindex = 1:numROIs
    if isempty(ROIdata.rois(ROIid(rindex)).rawdata) % replace whole vector
        ROIdata.rois(ROIid(rindex)).rawdata = Data(rindex, :);
        ROIdata.rois(ROIid(rindex)).rawneuropil = Neuropil(rindex, :);
    else
        ROIdata.rois(ROIid(rindex)).rawdata(FrameIndex) = Data(rindex, FrameIndex);
        ROIdata.rois(ROIid(rindex)).rawneuropil(FrameIndex) = Neuropil(rindex, FrameIndex);
    end
end
fprintf('\nFinished extracting signals from %d frames for %d ROIs in %.1f minutes\n', numFrames, numROIs, toc/60)

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
    fprintf('\tROIdata saved to: %s\n', saveFile);
end
