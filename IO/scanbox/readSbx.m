function [Images, Config, InfoFile] = readSbx(SbxFile, InfoFile, varargin)
% Loads 'Frames' of single .sbx file ('SbxFile'). Requires
% corresponding information file ('InfoFile').

LoadType = 'Direct'; % 'MemMap' or 'Direct'
Frames = 1:20; % indices of frames to load in 'Direct' mode, or 'all'
Channels = 1;
Depths = inf;
Verbose = true;

% Defaults
invert = true;
flipLR = false;
xavg = 4; %scanbox version 1 only

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
        
        % Determine frames to load
        if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
            Frames = 1:Config.Frames;
        elseif Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):Config.Frames];
        end
        numFramesToLoad = numel(Frames);
        Frames = Frames - 1; % offset b/c of "0" indexing
        
        % Determine channels to load
        if ischar(Channels) || (numel(Channels)==1 && Channels == inf)
            Channels = 1:Config.Channels;
        elseif Channels(end) == inf
            Channels = [Channels(1:end-2),Channels(end-1):Config.Channels];
        end
        numChannelsBeingRead = numel(Channels);
        
        % Determine depths to load
        if ischar(Depths) || (numel(Depths)==1 && Depths == inf)
            Depths = 1:Config.Depth;
        elseif Depths(end) == inf;
            Depths = [Depths(1:end-2),Depths(end-1):Config.Depth];
        end
        numDepths = numel(Depths);
        
        % Preallocate output
        if info.scanbox_version == 1
            Images = zeros(numChannelsBeingRead, Config.Width/xavg, Config.Height, numDepths, numFramesToLoad, 'uint16');
        elseif info.scanbox_version == 2
            Images = zeros(numChannelsBeingRead, Config.Width, Config.Height, numDepths, numFramesToLoad, 'uint16');
        end

        
        % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
        seekoperations = find(diff(Frames)~=1); %find any jumps within the frames to read out (jumps requiring seeking)
        if isempty(seekoperations) %no jumps
            numframesperread = numFramesToLoad; %all frames will be read in one read
            seekoperations = 1; %only one seek operation to first frame of FrameIndex
        else
            numframesperread = diff([0,seekoperations,numFramesToLoad]); %multiple reads required with various numbers of frames per read
            seekoperations = [1,seekoperations+1]; %indexes the first frame of each read within FrameIndex
        end
        
        % Open File
        info.fid = fopen(SbxFile);
        if(info.fid ~= -1)
            
            % Load Images
            if Verbose
                fprintf('Loading\t%d\tframe(s) from\t%s...', numel(Frames), SbxFile);
            end
            
            for index = 1:length(seekoperations)
                if(fseek(info.fid, Config.Height * Config.Width * Config.Depth * 2 * Frames(seekoperations(index)) * Config.Channels, 'bof')==0) % "2" b/c assumes uint16 => 2 bytes per record
                    temp = fread(info.fid, Config.Height * Config.Width * Config.Depth * Config.Channels * numframesperread(index), 'uint16=>uint16');
                    if info.scanbox_version == 1 % downsample/averaging
                        temp = mean(reshape(temp,[Config.Channels xavg Config.Width/xavg Config.Height Config.Depth numframesperread(index)]), 2);
                        temp = reshape(temp,[Config.numChannels Config.Width/xavg Config.Height Config.Depth numframesperread(index)]);
                    elseif info.scanbox_version == 2
                        temp = reshape(temp,[Config.Channels Config.Width Config.Height Config.Depth numframesperread(index)]);
                    end
                    Images(:,:,:,:,seekoperations(index):seekoperations(index)+numframesperread(index)-1) = temp(Channels,:,:,Depths,:); % save only requested channels & depths
                else
                    warning('fseek error...');
                    Images = [];
                end
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
                good = zeros(Config.Height, Config.Width, 1, numChannelsBeingRead, Config.Frames);
                for ii = 1:numChannelsBeingRead
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