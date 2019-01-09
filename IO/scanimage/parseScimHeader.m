function config = parseScimHeader(TifFile)
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

% config.version = 1;

%% Check input arguments
narginchk(0,1);
if ~exist('TifFile', 'var') || isempty(TifFile) % Prompt for file selection
    [TifFile, p] = uigetfile({'*.tif'}, 'Select ''tif'' files to load header info from:', directory);
    if isnumeric(TifFile)
        return
    end
    TifFile = fullfile(p, TifFile);
end


%% Set identifying info
config.type = 'scim';
config.FullFilename = TifFile;
[~, config.Filename, ~] = fileparts(TifFile);

%% Identify header information from file

% Load header
try
    header = scim_openTif(TifFile);
catch
    temp = imfinfo(TifFile);
    temp = parseSI5Header(temp(1).Software);
    header = temp.SI;
end
config.header = {header};

% Save important information
if isfield(header,'VERSION_MAJOR')
    config.Height = header.hRoiManager.linesPerFrame;
    config.Width = header.hRoiManager.pixelsPerLine;
    config.Depth = header.hStackManager.numSlices;
    config.ZStepSize = header.hStackManager.stackZStepSize;
    config.Channels = numel(header.hChannels.channelSave);
    config.FrameRate = header.hRoiManager.scanFrameRate;
    config.ZoomFactor = header.hRoiManager.scanZoomFactor;
    config.Frames = header.hStackManager.framesPerSlice;
elseif isfield(header,'scanimage') && header.scanimage.VERSION_MAJOR == 5
    config.Height = header.scanimage.linesPerFrame;
    config.Width = header.scanimage.pixelsPerLine;
    config.Depth = header.scanimage.stackNumSlices;
    config.ZStepSize = header.scanimage.stackZStepSize;
    config.Channels = size(header.scanimage.channelsSave,1);
    config.FrameRate = 1 / header.scanimage.scanFramePeriod;
    config.ZoomFactor = header.scanimage.zoomFactor;
    config.Frames = header.scanimage.acqNumFrames;
else
    config.Height = header.acq.linesPerFrame;
    config.Width = header.acq.pixelsPerLine;
    config.Depth = header.acq.numberOfZSlices;
    config.ZStepSize = header.acq.zStepSize;
    config.Channels = header.acq.numberOfChannelsSave;
    config.FrameRate = header.acq.frameRate;
    config.ZoomFactor = header.acq.zoomFactor;
    
    % Determine number of frames
    temp = imfinfo(TifFile);
    config.Frames = numel(temp)/(config.Channels*config.Depth);
end


%% DEFAULTS
% config.MotionCorrected = false;
% config.info = [];
config.Precision = 'uint16';
config.DimensionOrder = {'Height','Width','Channels','Frames'}; % default
config.Colors = {'green', 'red'};
config.size = sizeDimensions(config);


