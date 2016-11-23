function [Images, Config, InfoFile] = readSbx(SbxFile, InfoFile, varargin)
% Loads 'Frames' of single .sbx file ('SbxFile'). Requires
% corresponding information file ('InfoFile').

LoadType = 'Direct'; % 'MemMap' or 'Direct'
Frames = 1:20; % indices of frames to load in 'Direct' mode, or 'all'
Channels = 1;
Depths = inf;
Verbose = true;

% Defaults
invert = true;      % invert colormap boolean
FramesPerDepth = 1; % number of frames acquired at each depth before moving to next depth
flipLR = false;     % flip images across vertical axis
xavg = 4;           % number of pixels to average for each pixel (scanbox version 1 only)

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Type','type'}
                LoadType = varargin{index+1}; %'Direct' or 'MemMap'
                index = index + 2;
            case {'Frames','frames'} % indices of frames to load in 'Direct' mode
                Frames = varargin{index+1};
                index = index + 2;
            case {'Channels','channels'}
                Channels = varargin{index+1};
                index = index + 2;
            case {'Depths', 'depths', 'Depth', 'depth'}
                Depths = varargin{index+1};
                index = index + 2;
            case {'Invert', 'invert'}
                invert = varargin{index+1};
                index = index + 2;
            case {'Flip', 'flip'}
                flipLR = varargin{index+1};
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

if ~exist('SbxFile', 'var') || isempty(SbxFile)
    directory = CanalSettings('DataDirectory');
    [f,p] = uigetfile({'*.sbx'}, 'Choose scanbox file to load', directory);
    if isnumeric(f)
        Images = []; return
    end
    SbxFile = fullfile(p,f);
end

if ~exist('InfoFile', 'var') || isempty(InfoFile)
    InfoFile = sbxIdentifyFiles(SbxFile);
    InfoFile = InfoFile{1};
end

%% Load In Acquisition Information
load(InfoFile, 'info'); %load in 'info' variable

% Parse acquisition information
Config = parseSbxHeader(InfoFile);
Config.loadType = LoadType;

%% Load In Images
switch LoadType
    
    case {'MemMap', 'memmap'}
        
        Images = MappedTensor(SbxFile, Config.size, 'Class', Config.Precision);
        Images = permute(Images, [3 2 5 1 4]);
        % Images = intmax('uint16') - Images; % computation doesn't carry over
        Config.DimensionOrder = Config.DimensionOrder([3 2 5 1 4]);
        Config.size = size(Images);
        % Images = memmapfile(SbxFile,...
        %     'Format', {'uint16', [info.numChannels info.Width info.Height info.numFrames], 'Frames'}); %format data as best as possible into individual frames
        % 'Offset', nsamples*2*info.nchan,... %skip first frame as it's incomplete
        % Frames = permute(intmax('uint16')-x.Data.Frames,[3,2,5,1,4]);
        
    case {'Direct', 'direct'}
        
        % Determine encoding
        switch Config.Precision
            case 'uint16'
                BytesPerPixel = 2;
        end
        PixelsPerFrame = Config.Height*Config.Width;
        
        % Determine frames to load
        if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
            Frames = 1:Config.Frames;
        elseif Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):Config.Frames];
        end
        numFrames = numel(Frames);
        
        % Determine channels to load
        if ischar(Channels) || (numel(Channels)==1 && Channels == inf)
            Channels = 1:Config.Channels;
        elseif Channels(end) == inf
            Channels = [Channels(1:end-2),Channels(end-1):Config.Channels];
        end
        numChannels = numel(Channels);
        
        % Determine depths to load
        if ischar(Depths) || (numel(Depths)==1 && Depths == inf)
            Depths = 1:Config.Depth;
        elseif Depths(end) == inf;
            Depths = [Depths(1:end-2),Depths(end-1):Config.Depth];
        end
        numDepths = numel(Depths);
        
        % Determine indices within file to load
        PacketSize = FramesPerDepth * numDepths;            % number of frames in cycle
        F = reshape(1:PacketSize,FramesPerDepth,numDepths); % frame index corresponding to frame 1:FramesPerDepth of each depth
        Findex = rem(Frames-1,FramesPerDepth)+1;            % index of each requested frame within FramesPerDepth
        Frames = bsxfun(@plus, floor((Frames-1)/FramesPerDepth)'*PacketSize,F(Findex,:))'; % frame ID within whole movie
        [Frames,order] = sort(Frames(:));                   % list of frame indices to load
          
        % Preallocate output
        if info.scanbox_version ~= 1
            Images = zeros(numChannels*PixelsPerFrame*FramesPerDepth*numDepths*floor(numFrames/FramesPerDepth),1, 'uint16');
        else % version 1
            Images = zeros(numChannels,Config.Width/xavg,Config.Height,numDepths,numFrames, 'uint16');
        end

        % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
        jumps = diff([0;Frames])-1;     % number of frames to jump
        startID = find(jumps);          % indices of jump locations relative to frame indices
        if ~jumps(1)                    % loading very first frame in file
            startID = [1;startID];      % ensures chunk containing first frame is loaded
        end
        jumps = jumps(startID);
        N = diff([startID-1;numFrames*numDepths]); % number of frames to load in on each loop

        
        % Open File
        info.fid = fopen(SbxFile);
        if(info.fid ~= -1)
            
            % Load Images
            if Verbose
                fprintf('Loading\t%d\tframe(s) from\t%s...', numFrames, SbxFile);
            end
            
            % Cycle through loading in frames
            frewind(info.fid); % reset file
            for index = 1:length(startID)
                if(fseek(info.fid, PixelsPerFrame*BytesPerPixel*Config.Channels*jumps(index),'cof')==0)
                    if info.scanbox_version ~= 1
                        Images((startID(index)-1)*PixelsPerFrame+1:PixelsPerFrame*(startID(index)+N(index)-1)) = fread(info.fid, PixelsPerFrame*Config.Channels*N(index), 'uint16=>uint16');
                    else % downsample/averaging
                        temp = fread(info.fid, PixelsPerFrame * Config.Depth * Config.Channels * N(index), 'uint16=>uint16');
                        temp = mean(reshape(temp,[Config.Channels xavg Config.Width/xavg Config.Height Config.Depth N(index)]), 2);
                        temp = reshape(temp,[Config.numChannels Config.Width/xavg Config.Height Config.Depth N(index)]);
                        Images(:,:,:,:,startID(index):startID(index)+N(index)-1) = temp(Channels,:,:,Depths,:); % save only requested channels & depths
                    end
                else
                    warning('fseek error...');
                    Images = [];
                end
            end
            
            % Reshape images & keep only desired channels
            if info.scanbox_version ~= 1
                if FramesPerDepth == 1 || numDepths == 1
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,numDepths,numFrames); % reshape vector to usable format
                else % FramesPerDepth > 1 -> requested frames could be out of order
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,numel(Frames));
                    warning('Never tested that the next line works, please verify this is valid...');
                    Images = Images(:,:,:,order);
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,numDepths,numFrames); % reshape vector to usable format
                end
                Images = Images(Channels,:,:,:,:); % keep only requested channels
            end
            
            % Update dimensions (v.1 only)
            if info.scanbox_version == 1
                Config.Width = Config.Width/xavg;
            end
            
            % Reorder dimensions
            Images = permute(Images, [3 2 4 1 5]); % flip colormap and reorder (original [1,3,2,4] => rotate images)
            
            % Correct for nonuniform spatial sampling
            if info.scanbox_version == 1
                S = sparseint;
                info.Width = size(S,2);
                good = zeros(Config.Height, Config.Width, 1, numChannels, Config.Frames);
                for ii = 1:numChannels
                    for iii = 1:1
                        for iiii = 1:Config.Frames
                            good(:,:,iii,ii,iiii) = Images(:,:,iii,ii,iiii) * S; % correct for non-uniform sampling
                        end
                    end
                end
                Images = good;
            end

            % Flip colormap
            if invert
                Images = intmax('uint16') - Images;
            end
            
            if flipLR
                Images = fliplr(Images);
            end
            
        else
            warning('unable to open file: %s', SbxFile);
            Images = [];
        end
        fclose(info.fid);
        
        Config.DimensionOrder = Config.DimensionOrder([3 2 4 1 5]);
        
        Config.size = size(Images);
        
        if Verbose
            fprintf('\tComplete\n');
        end

end