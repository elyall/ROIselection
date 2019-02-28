function Images = readTiff(TifFile, varargin)
% Loads 'Frames' of single .tif file ('TifFile')

Frames = [1,inf]; % indices of frames to load in 'Direct' mode, or 'all'
AutoScale = false;
Verbose = false;

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Frames','frames'} % indices of frames to load in 'Direct' mode
                Frames = varargin{index+1};
                index = index + 2;
            case {'AutoScale','Scale'}
                AutoScale = varargin{index+1};
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
    [TifFile,p] = uigetfile({'*.tif'}, 'Choose scanbox file to load');
    if isnumeric(TifFile)
        Images = []; return
    end
    TifFile = fullfile(p,TifFile);
end


%% Load in header
warning('off', 'MATLAB:imagesci:tifftagsread:nextIfdPointerOutOfRange');
info = imfinfo(TifFile,'tif');
totalFrames = numel(info);


%% Determine frames to load
if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
    Frames = 1:totalFrames;
elseif Frames(end) == inf
    Frames = [Frames(1:end-2),Frames(end-1):totalFrames];
end
numFrames = numel(Frames);


%% Determine data type
Class = 'double';
switch info(1).SampleFormat
    case 'Unsigned integer'
        if info(1).BitDepth==8
            Class = 'uint8';
        elseif info(1).BitDepth==16
            Class = 'uint16';
        elseif info(1).BitDepth==32
            Class = 'uint32';
        end
end
    

%% Load in frames
if Verbose
    fprintf('Loading\t%d\tframe(s) from\t%s...', numFrames, TifFile);
    Fn=parfor_progress(numFrames);
end

% Initialize output
switch info(1).ColorType
    case 'truecolor'
        Images = zeros(info(1).Height, info(1).Width, 3, numFrames, Class);
    otherwise
        Images = zeros(info(1).Height, info(1).Width, 1, numFrames, Class);
end

% Read in file
tif=Tiff(TifFile,'r');
for findex = 1:numFrames
    tif.setDirectory(Frames(findex))
    Images(:,:,:,findex)=tif.read();
    if Verbose
        parfor_progress(Fn);
    end
end
tif.close;

if Verbose
    parfor_progress(Fn,0);
end

if AutoScale
    Images = double(Images);
    Images = Images*double(intmax(Class))/max(Images(:));
    Images = cast(Images,Class);
end
