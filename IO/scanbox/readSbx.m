function [Images, Config, InfoFile] = readSbx(SbxFile, varargin)
%READSBX   Load .sbx files
%   Images = readSbx() will prompt user to select a file to load.
%
%   [Images, Config, InfoFile] = readSbx(SbxFile) loads in the first
%   default number of frames from the sbx file input, and returns the
%   metadata and filename of associated info file.
%
%   Images = readSbx(..., 'Frames', X) loads in specific frames from file.
%   Set X equal to inf to load in all frames.
%
%   Images = readSbx(..., 'Depths', Y) loads in specific depths from file.
%   Set Y equal to inf to load in all depths.
%
%   Images = readSbx(..., 'Channels', Z) loads in all Channels, but then
%   throws out all Channels except those requested (would take much much
%   longer to just load a specific channel, as channels are interleaved
%   every 16 bits).


% Default parameters that can be adjusted
LoadType = 'Direct';% specifies how the data is loaded, can be: 'MemMap' or 'Direct'
Frames = 1:20;      % indices of frames to load in 'Direct' mode, or 'all'
Channels = 1;       % default channels to load
Depths = inf;       % default depths to load
Verbose = true;     % booleon determining whether to display progress bar
invert = true;      % invert colormap boolean
FramesPerDepth = 1; % number of frames acquired at each depth before moving to next depth
flipLR = false;     % flip images across vertical axis
xavg = 4;           % scanbox version 1: number of pixels to average for each pixel

% Placeholders
InfoFile = '';


%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'InfoFile'
                InfoFile = varargin{index+1};
                index = index + 2;
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

if isempty(InfoFile)
    InfoFile = sbxIdentifyFiles(SbxFile);
    InfoFile = InfoFile{1};
end

if iscolumn(Frames)
    Frames = Frames';
end


%% Load In Acquisition Information
load(InfoFile, 'info'); %load in 'info' variable

% Parse acquisition information
Config = parseSbxHeader(InfoFile);
if info.scanbox_version == 1 && Config.Depth > 1
    error('Never coded proper loading of multiple depths for old data');
end
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
        HxW = Config.Height*Config.Width;
        
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
        if Config.Depth>1 % if Config.Depth==1 then frame indices is the same as the file indices
            PacketSize = FramesPerDepth * Config.Depth;             % number of frames in cycle
            F = reshape(1:PacketSize,FramesPerDepth,Config.Depth);	% frame index corresponding to frame 1:FramesPerDepth of each depth
            Findex = rem(Frames-1,FramesPerDepth)+1;                % index of each requested frame within FramesPerDepth
            Frames = bsxfun(@plus, floor((Frames-1)/FramesPerDepth)'*PacketSize,F(Findex,:)); % frame ID within whole movie
            Frames = Frames(:,Depths)';                             % load only requested depths; transpose for vectorizing
            [Frames,order] = sort(Frames(:));                       % list of frame indices to load
        else
            Frames = Frames';
        end
          
        % Preallocate output
        if info.scanbox_version ~= 1
            Images = zeros(numChannels*HxW*FramesPerDepth*numDepths*numFrames,1, 'uint16');
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
                fprintf('Loading %7.d frame(s) & %2.d depth(s) from %s...', numFrames, numDepths, SbxFile);
                t=tic;
            end
            
            % Cycle through loading in frames
            frewind(info.fid); % reset file
            for index = 1:length(startID)
                if(fseek(info.fid, HxW*BytesPerPixel*Config.Channels*jumps(index),'cof')==0)
                    if info.scanbox_version ~= 1
                        Images((startID(index)-1)*HxW+1:HxW*(startID(index)-1+N(index))) = fread(info.fid, HxW*Config.Channels*N(index), 'uint16=>uint16');
                    else % downsample/averaging to avoid overloading memory
                        temp = fread(info.fid, HxW*Config.Channels*N(index), 'uint16=>uint16');
                        temp = mean(reshape(temp,[Config.Channels xavg Config.Width/xavg Config.Height 1 N(index)]), 2);
                        temp = reshape(temp,[Config.numChannels Config.Width/xavg Config.Height 1 N(index)]);
                        Images(:,:,:,:,startID(index):startID(index)+N(index)-1) = temp(Channels,:,:,1,:); % save only requested channels & depths
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
                    warning('Never tested that the next line works, please verify it is valid...');
                    Images = Images(:,:,:,order); % <= VERIFY THIS LINE REORDERS DATA CORRECTLY (really only necessary when rem(numFrames,FramesPerDepth)~=0)
                    Images = reshape(Images,numChannels,Config.Width,Config.Height,numDepths,numFrames); % reshape vector to usable format
                end
                Images = Images(Channels,:,:,:,:); % keep only requested channels
            end
            
            % Reorder dimensions
            Images = permute(Images, [3 2 4 1 5]); % flip and move channels to fourth dimension
            
            % Correct for nonuniform spatial sampling
            if info.scanbox_version == 1
                Config.Width = Config.Width/xavg; % update dimensions for Config output
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
            
            % Flip across Y axis
            if flipLR
                Images = fliplr(Images);
            end
            
        else
            warning('unable to open file: %s', SbxFile);
            Images = [];
        end
        fclose(info.fid); % close file
        
        
        Config.DimensionOrder = Config.DimensionOrder([3 2 4 1 5]); % update dimensions to that of output as compared to what they are in the file
        Config.size = size(Images); % update size of output
        
        
        if Verbose
            t=toc(t);
            if t<60
                s = 'seconds';
            elseif t<3600
                t = t/60;
                s = 'minutes';
            else
                t = t/3600;
                s = 'hours -> probably overloading RAM and using swap space, it''s recommended to load less frames at a time';
            end
            fprintf('\tComplete (%.2f %s)\n',t,s);
        end

end