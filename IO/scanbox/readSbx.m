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
        PacketSize = FramesPerDepth * numDepths;
        F = reshape(1:PacketSize,FramesPerDepth,numDepths);
        Findex = rem(Frames-1,FramesPerDepth)+1;
        Frames = bsxfun(@plus, floor(Frames/numDepths)',F(Findex,:)); % STILL WRONG for 
        
        
        Frames = Frames-1; % offset b/c of "0" indexing when loading

        
        
        if rem(numFrames,FramesPerDepth)~=0
            Frames = [Frames,nan(1,FramesPerDepth-rem(numFrames,FramesPerDepth))];
        end
        Frames = reshape(Frames, FramesPerDepth, ceil(numFrames/FramesPerDepth));
        Frames = repmat(Frames, [1,1,FramesPerDepth]);
        Frames = bsxfun(@plus, repmat(temp',1,1,FramesPerDepth)*Config.Depth, 1:FramesPerDepth:Config.Depth*FramesPerDepth);
        Frames = bsxfun(@plus, Frames, permute(0:FramesPerDepth-1,[1,3,2]));
        
        % Preallocate output
        if info.scanbox_version ~= 1
            Images = zeros(numChannels*PixelsPerFrame*FramesPerDepth*numDepths*floor(numFrames/FramesPerDepth), 'uint16');
        else % version 1
            Images = zeros(numChannels,Config.Width/xavg,Config.Height,numDepths,numFrames, 'uint16');
        end

        % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
        seekoperations = find(diff(Frames)~=1); %find any jumps within the frames to read out (jumps requiring seeking)
        if isempty(seekoperations) %no jumps
            numframesperread = numFrames; %all frames will be read in one read
            seekoperations = 1; %only one seek operation to first frame of FrameIndex
        else
            numframesperread = diff([0,seekoperations,numFrames]); %multiple reads required with various numbers of frames per read
            seekoperations = [1,seekoperations+1]; %indexes the first frame of each read within FrameIndex
        end
        
        % Open File
        info.fid = fopen(SbxFile);
        if(info.fid ~= -1)
            
            % Load Images
            if Verbose
                fprintf('Loading\t%d\tframe(s) from\t%s...', numFrames, SbxFile);
            end
            
            % Cycle through loading in frames
            for index = 1:length(seekoperations)
                if(fseek(info.fid, PixelsPerFrame * BytesPerPixel * Config.Depth * Config.Channels * Frames(seekoperations(index)), 'bof')==0)
                    if info.scanbox_version ~= 1
                        Images(seekoperations(index):seekoperations(index)+numframesperread(index)-1) = fread(info.fid, PixelsPerFrame * Config.Depth * Config.Channels * numframesperread(index), 'uint16=>uint16');
                    else % downsample/averaging
                        temp = fread(info.fid, PixelsPerFrame * Config.Depth * Config.Channels * numframesperread(index), 'uint16=>uint16');
                        temp = mean(reshape(temp,[Config.Channels xavg Config.Width/xavg Config.Height Config.Depth numframesperread(index)]), 2);
                        temp = reshape(temp,[Config.numChannels Config.Width/xavg Config.Height Config.Depth numframesperread(index)]);
                        Images(:,:,:,:,seekoperations(index):seekoperations(index)+numframesperread(index)-1) = temp(Channels,:,:,Depths,:); % save only requested channels & depths
                    end
                else
                    warning('fseek error...');
                    Images = [];
                end
            end
            
            % Reshape images & keep only desired channels
            if info.scanbox_version ~= 1
                if FramesPerDepth == 1
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,numDepths,numFrames); % reshape vector to usable format
                else
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,FramesPerDepth,numDepths,floor(numFrames/FramesPerDepth)); % reshape vector to usable format
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