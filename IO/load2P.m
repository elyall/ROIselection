function [Images, loadObj, Config] = load2P(ImageFiles, varargin)
%LOAD2P    Loads multiple .sbx or .tif files.
%   IMAGES = load2P() will prompt user to select a sbx file to load and
%   returns the first 20 frames from that file.
%
%   [IMAGES, LOADOBJ, CONFIG] = load2P(IMAGEFILES) loads in the first 20
%   frames from IMAGEFILES, the location of a single file or a cell array
%   of strings containing multiple filenames. LOADOBJ is a struct detailing
%   the data loaded in IMAGES. CONFIG is the metadata for all files located
%   in IMAGEFILES.
%
%   [...] = load2P(..., 'Frames', FRAMES) a vector that specifies the exact
%   frame indices to load from the file. Set the end of X equal to inf to
%   load to the end of the file. (default = 1:20)
%
%   [...] = load2P(..., 'Depths', DEPTHS) a vector that specifies the
%   depths to load from the file. (default = inf)
%
%   [...] = load2P(..., 'Channels', CHANNELS) a vector that specifies which
%   channels to include in the output. (default = 1) (Note: all channels
%   are initially loaded into memory due to the channels being interleaved
%   at the 16-bit scale causing a burden to do so otherwise)
%
%   [...] = load2P(..., 'IndexType', TYPE) 'absolute' or 'relative' that
%   sets whether the frame indices specified are absolute indices, or
%   relative to the depths requested. (default = 'absolute')
%
%   [...] = load2P(..., 'Save', X) empty string to not save images to a
%   file, true to prompt for a filename to save to, or a string specifying
%   the filename to save the images to (will also prompt user if file
%   already exists). (default = '')
%
%   [...] = load2P(..., 'Type', LOADTYPE) 'Direct' or 'MemMap' specifying
%   whether the frames should be memory mapped or loaded directly into
%   memory. (default = 'Direct')
%
%   [...] = load2P(..., 'double') forces IMAGES to be of class 'double'.
%
%   [...] = load2P(..., 'verbose') displays loading bar.
%

% Default parameters that can be adjusted
LoadType = 'Direct';    % 'MemMap' or 'Direct' 
Frames = 1:20;          % indices of frames to load in 'Direct' mode
IndexType = 'absolute';	% 'absolute' or 'relative' -> specifies whether 'Frames' indexes the absolute frame index or the relative frame index for the depth(s) requested (doesn't matter if only 1 depth in the file)
Channels = 1;           % default channels to load
Depths = inf;           % default depths to load
Double = false;         % booleon determining whether output images are forced to be of class double
Verbose = false;        % booleon determining whether to display progress bar
SaveFile = '';          % empty string, True to prompt for filename selection, or filename of file to save images to

% Placeholders
directory = cd; % default directory when prompting user to select a file

%% Initialize Parameters
if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose images file(s) to load', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles)
        Images = []; return
    elseif iscellstr(ImageFiles)
        for index = 1:numel(ImageFiles)
            ImageFiles{index} = fullfile(p,ImageFiles{index});
        end
    else
        ImageFiles = {fullfile(p,ImageFiles)};
    end
elseif ischar(ImageFiles)
    ImageFiles = {ImageFiles};
elseif isstruct(ImageFiles)
    loadObj = ImageFiles;
end

if iscellstr(ImageFiles) && any(cellfun(@isempty, regexp(ImageFiles,'.[/w.sbx/w.tif]$')))
    error('Make sure all files input are either .sbx or .tif files');
end
numFiles = numel(ImageFiles);

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case isstruct(varargin{index})
                loadObj = varargin{index};
                index = index + 1;
            case {'Type','type'}
                LoadType = varargin{index+1};
                index = index + 2;
            case {'Frames','frames','Frame','frame'} % indices of frames to load in 'Direct' mode
                Frames = varargin{index+1};
                index = index + 2;
            case 'IndexType'
                IndexType = varargin{index+1};
                index = index + 2;
            case {'Channels', 'channels', 'Channel', 'channel'}
                Channels = varargin{index+1};
                index = index + 2;
            case {'Depths', 'depths', 'Depth', 'depth'}
                Depths = varargin{index+1};
                index = index + 2;
            case {'Save','save','saveFile'}
                SaveFile = varargin{index+1};
                index = index + 2;
            case {'Double', 'double'}
                Double = ~Double;
                index = index + 1;
            case {'Verbose', 'verbose'}
                Verbose = ~Verbose;
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

% Check LoadType
if ~any(strcmp(LoadType,{'Direct','MemMap'}))
    error('Load type ''%s'' not recognized. Load type has to be either ''Direct'' or ''MemMap''', LoadType);
end

% Update frames to load
if ischar(Frames) && strcmp(Frames, 'all')
    Frames = [1 inf];
elseif isnumeric(Frames) && isvector(Frames) && size(Frames,1)~=1
    Frames = Frames';
end

if ~exist('loadObj', 'var')
    % Initialize LoadObj
    for findex = 1:numFiles
        loadObj.files(findex).FullFilename = ImageFiles{findex};
        loadObj.files(findex).Frames = [];
        loadObj.files(findex).ext = '';
    end
else
    LoadType = loadObj.Type; % set load type based on loadObj
end

% Update LoadObj
loadObj.Type = LoadType;
loadObj.FrameIndex = [];


%% Load in header information
Config = load2PConfig(ImageFiles);
for findex = 1:numFiles
    loadObj.files(findex).Frames = Config(findex).Frames;
end


%% Load images
switch LoadType
    
%% Load Direct
    case 'Direct'

        % Determine Dimensions are Equal Across Files
        if range([Config(:).Height]) ~= 0 || range([Config(:).Width]) ~= 0 || range([Config(:).Depth]) ~= 0
            error('Data need to be the same size...');
        end
        % Scanbox version 1 fix
        if strcmp(Config(1).type,'sbx') && Config(1).header{1}.scanbox_version == 1
            Config.Width = 796;
        end
        
        % Determine number of frames in each file and which frames to load
        % from each file
        numFrames = [Config(:).Frames];
        FrameIndex = cell(numFiles, 1);
        if iscell(Frames) % specific frames from each file are designated
            if numel(Frames) == 1 && numFiles > 1 % single cell entry specified
                Frames = repmat(Frames, numFiles, 1);
            elseif numel(Frames) ~= numFiles
                error('Must specify which frames to load from each file');
            end
            for findex = 1:numFiles
                FrameIndex{findex} = Frames{findex};
                if FrameIndex{findex}(end) == inf
                    try
                        FrameIndex{findex} = cat(2, FrameIndex{findex}(1:end-1), FrameIndex{findex}(end-1)+1:numFrames(findex));
                    catch % only 'inf' input (numel(Frames)==1)
                        FrameIndex{findex} = 1:numFrames(findex);
                    end
                end
            end
        elseif isnumeric(Frames) % take all files to be one experiment and select relative frames
            cumsumFrames = [0,cumsum(numFrames)];
            totalFrames = cumsumFrames(end);
            cumsumFrames(end) = [];
            if Frames(end) == inf
                try
                    Frames = cat(2, Frames(1:end-1), Frames(end-1)+1:totalFrames);
                catch % only 'inf' input (numel(Frames)==1)
                    Frames = 1:totalFrames;
                end
            end
            temp = bsxfun(@minus, Frames, cumsumFrames');
            for findex = 1:numFiles
                FrameIndex{findex} = temp(findex, temp(findex,:)>=1 & temp(findex,:)<=numFrames(findex));
            end
        end
        
        % Determine Channels to Load (FILES ALL MUST HAVE DESIRED CHANNEL
        % OR WILL ERROR)
        if ischar(Channels) || (numel(Channels)==1 && Channels == inf)
            Channels = 1:min([Config(:).Channels]); % load all channels (based on file with minimum # of channels)
        elseif Channels(end) == inf
            Channels = [Channels(1:end-2),Channels(end-1):min([Config(:).Channels])];
        end
        numChannels = numel(Channels);
        
        % Determine Depths to load (defaults to minimum # of channels)
        if ischar(Depths) || (numel(Depths)==1 && Depths == inf)
            Depths = 1:min([Config(:).Depth]); % load all channels (based on file with minimum # of depths)
        elseif Depths(end) == inf
            Depths = [Depths(1:end-2),Depths(end-1):min([Config(:).Depths])];
        end
        numDepths = numel(Depths);
        
        % Determine absolute frame indices to load from each file
        numFrames = zeros(1,numFiles);
        for index = 1:numFiles
            if  ~isempty(FrameIndex{index})
                if Config(index).Depth>1
                    [depthID,relativeIndex] = idDepth(Config(index),[],'IndexType',IndexType,'Frames',FrameIndex{index},'Depths',Depths,'FramesPerDepth',Config(index).FramesPerDepth); % determine file indices of frames requested
                    numFrames(index) = numel(relativeIndex);
                    numDepths = sum(any(~isnan(depthID),1));
                    depthID = depthID';
                    FrameIndex{index} = sort(depthID(:))'; % list of frame indices to load
                    FrameIndex{index}(isnan(FrameIndex{index})) = []; % remove NaN's
                    loadObj.FrameIndex = cat(1, loadObj.FrameIndex, cat(2, index*ones(numFrames(index),1), relativeIndex, depthID'));
                else
                    numDepths = 1;
                    numFrames(index) = numel(FrameIndex{index});
                    loadObj.FrameIndex = cat(1, loadObj.FrameIndex, cat(2, index*ones(numel(FrameIndex{index}),1), FrameIndex{index}'));
                end
            end
        end
        
        % Load Images
        Images = zeros(Config(1).Height, Config(1).Width, numDepths, numChannels, sum(numFrames), 'uint16');
        startFrame = cumsum([1,numFrames(1:end-1)]);
        for index = 1:numFiles
            [~,~,loadObj.files(index).ext] = fileparts(ImageFiles{index});
            if ~isempty(FrameIndex{index})
                switch loadObj.files(index).ext
                    case '.sbx'
                        Images(:,:,:,:,startFrame(index):startFrame(index)+numFrames(index)-1)...
                            = readSbx(ImageFiles{index}, 'Type', 'Direct', 'Frames', FrameIndex{index}, 'IndexType', 'absolute', 'Channels', Channels, 'Depth', Depths, 'FramesPerDepth', Config(index).FramesPerDepth, 'Verbose', Verbose);
                    case '.tif'
                        Images(:,:,:,:,startFrame(index):startFrame(index)+numFrames(index)-1)...
                            = readScim(ImageFiles{index}, 'Frames', FrameIndex{index}, 'Channels', Channels, 'Verbose', Verbose);
                    case '.imgs'
                        Images(:,:,:,:,startFrame(index):startFrame(index)+numFrames(index)-1)...
                            = readImgs(ImageFiles{index}, 'Type', 'Direct', 'Frames', FrameIndex{index}, 'Channels', Channels);
                end
            end
            
        end
        
        if Double && ~isa(Images, 'double')
            Images = double(Images);
        end

        
%% Load MemMap
    case 'MemMap'
        
        if numFiles > 1
            warning('Cannot load more than one file with MemMap. Loading first file...');
            % numFiles = 1;
            Config = Config(1);
            ImageFiles = ImageFiles(1);
            loadObj.files(2:end) = [];
        end
        
        [~,~,loadObj.files.ext] = fileparts(ImageFiles{1});
        switch loadObj.files.ext
            case '.sbx'
                Images = readSbx(ImageFiles{1}, [], 'Type', LoadType);
            case '.tif'
                Images = readScim(ImageFiles, 'Type', LoadType);
            case '.imgs'
                loadObj.files.memmap = readImgs(ImageFiles, 'Type', LoadType);
                Images = loadObj.files.memmap.Data.Frames;
        end
        loadObj.FrameIndex = cat(2, ones(Config.Frames, 1), (1:Config.Frames)');
end


%% Update data information
loadObj.Precision = class(Images(1));
[loadObj.Height, loadObj.Width, loadObj.Depth, loadObj.Channels, loadObj.Frames] = size(Images);
loadObj.DimensionOrder = {'Height', 'Width', 'Depth', 'Channels', 'Frames'};
loadObj.size = [loadObj.Height, loadObj.Width, loadObj.Depth, loadObj.Channels, loadObj.Frames];

loadObj.FrameRate = mode([Config(:).FrameRate]);


%% Save Images to file
if ~isempty(SaveFile)
    
    % Determine filename
    if isequal(SaveFile,true) || exist('SaveFile', 'file') % if file previously exists, prompt for filename
        if isequal(SaveFile,true) % set filename to be default filename
            [p,f,~] = fileparts(ImageFiles{1});
            SaveFile = fullfile(p, [f,'.mat']); 
        end
        [SaveFile, p] = uiputfile({'*.mat;*.tif;*.sbx'}, 'Save images as:', SaveFile); % have user select the file
        if isnumeric(SaveFile)  % user hit cancel
            return              % return outputs without saving images to file
        end
        SaveFile = fullfile(p,SaveFile);
    end
    
    % Save images to file
    [~,~,ext] = fileparts(SaveFile);
    switch ext
        case '.mat'
            save(SaveFile, 'Images', 'Config', 'numFrames', 'info', '-v7.3');
        otherwise
            save2P(SaveFile,Images);
    end
    
end

