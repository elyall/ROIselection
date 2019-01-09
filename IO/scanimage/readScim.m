function [Images, Config] = readScim(TifFile, varargin)
% Loads 'Frames' of single .tif file ('TifFile')

Frames = inf; % indices of frames to load in 'Direct' mode, or 'all'
Channels = 1;
Depths = 1;
Verbose = true;

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Frames','frames'} % indices of frames to load in 'Direct' mode
                Frames = varargin{index+1};
                index = index + 2;
            case {'Channels','channels'}
                Channels = varargin{index+1};
                index = index + 2;
            case {'Depth','depth','Depths','depths'}
                Depths = varargin{index+1};
                index = index + 2;
            case {'Verbose', 'verbose'}
                Verbose = varargin{index+1};
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

if ~exist('TifFile', 'var') || isempty(TifFile)
    [TifFile,p] = uigetfile({'*.sbx'}, 'Choose scanbox file to load');
    if isnumeric(TifFile)
        Images = []; return
    end
    TifFile = fullfile(p,TifFile);
end

%% Load In Acquisition Information
Config = parseScimHeader(TifFile);

%% Load In Images

% Determine channels to load
if ischar(Channels) || (numel(Channels)==1 && Channels == inf)
    Channels = 1:Config.Channels;
elseif Channels(end) == inf
    Channels = [Channels(1:end-2),Channels(end-1):Config.Channels];
end

% Determine depths to load
if ischar(Depths) || (numel(Depths)==1 && Depths == inf)
    Depths = 1:Config.Depth;
elseif Depths(end) == inf
    Depths = [Depths(1:end-2),Depths(end-1):Config.Depth];
end

% Determine frames to load
if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
    Frames = 1:Config.Frames;
elseif Frames(end) == inf
    Frames = [Frames(1:end-2),Frames(end-1):Config.Frames];
end

% Determine frame indices to load
info = imfinfo(TifFile,'tif');
totalFrames = numel(info);
IDs = idDepth(Config.Depth*Config.Channels,totalFrames);
IDs = reshape(IDs,[size(IDs,1),Config.Channels,Config.Depth]);
IDs = IDs(Frames,Channels,Depths);
IDs = permute(IDs,[2,3,1]);
IDs = IDs(:);

% Load in images
if Verbose
    fprintf('Loading\t%d\tframe(s) from\t%s...', numel(IDs), TifFile);
end
Images = readTiff(TifFile, 'Frames', IDs);

% Rearrange to correct dimensions
dim = size(Images);
Images = reshape(Images,[dim(1:2),numel(Channels),numel(Depths),numel(Frames)]);
Images = permute(Images, [1,2,4,3,5]);

if Verbose
    fprintf('\tComplete\n');
end
        

