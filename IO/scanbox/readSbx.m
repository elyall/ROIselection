function [Images, Config, InfoFile] = readSbx(SbxFile, varargin)
%READSBX   Loads .sbx files
%   IMAGES = readSbx() will prompt user to select a sbx file to load and
%   returns the first 20 frames from that file.
%
%   [IMAGES, CONFIG, INFOFILE] = readSbx(SBXFILE) loads in the first 20
%   frames from the sbx file SBXFILE, and returns the metadata and
%   filename of the metadata file.
%
%   [...] = readSbx(..., 'Frames', FRAMES) a vector that specifies the
%   exact frame indices to load from the file. Set the end of X equal to
%   inf to load to the end of the file. (default = 1:20)
%
%   [...] = readSbx(..., 'Depths', DEPTHS) a vector that specifies the
%   depths to load from the file. (default = inf)
%
%   [...] = readSbx(..., 'Channels', CHANNELS) a vector that specifies
%   which channels to include in the output. (default = 1) (Note: all
%   channels are initially loaded into memory due to the channels being
%   interleaved at the 16-bit scale causing a burden to do so otherwise)
%
%   [...] = readSbx(..., 'IndexType', TYPE) 'absolute' or 'relative' that
%   sets whether the frame indices specified are absolute indices, or
%   relative to the depths requested. (default = 'absolute')
%
%   [...] = readSbx(..., 'InfoFile', INFOFILENAME) specifies a specific
%   metadata file to use. (default = [SBXFILE(1:end-3),'mat'])
%
%   [...] = readSbx(..., 'Type', LOADTYPE) 'Direct' or 'MemMap' specifying
%   whether the frames should be memory mapped or loaded directly into
%   memory. (default = 'Direct')
%
%   [...] = readSbx(..., 'fliplr') flips images across the vertical axis.
%
%   [...] = readSbx(..., 'invert') inverts the colormap.
%
%   [...] = readSbx(..., 'verbose') displays loading bar. (default=on)
%
%   [...] = readSbx(..., 'organizeDepths') de-interleaves different depths.
%   (default=on)
%


% Default parameters that can be adjusted
LoadType = 'Direct';    % specifies how the data is loaded, can be: 'MemMap' or 'Direct'
Frames = 1:20;          % indices of frames to load in 'Direct' mode
IndexType = 'absolute'; % 'absolute' or 'relative' -> specifies whether 'Frames' indexes the absolute frame index or the relative frame index for the depth(s) requested (doesn't matter if only 1 depth in the file)
Channels = 1;           % default channels to load
Depths = inf;           % default depths to load
FramesPerDepth = 1;     % specifies number of frames taken at given depth before moving on to next depth
verbose = true;         % booleon determining whether to display progress bar
invert = false;         % invert colormap boolean
flipLR = false;         % flip images across vertical axis
organizeDepths=true;    % booleon determining whether to reshape file with multiple depths into 5D matrix or leave frames interleaved
xavg = 4;               % scanbox version 1: number of pixels to average for each pixel

% Placeholders
InfoFile = '';
directory = cd;

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
            case 'IndexType'
                IndexType = varargin{index+1};
                index = index + 2;
            case {'Channels','channels'}
                Channels = varargin{index+1};
                index = index + 2;
            case {'Depths', 'depths', 'Depth', 'depth'}
                Depths = varargin{index+1};
                index = index + 2;
            case 'FramesPerDepth'
                FramesPerDepth = varargin{index+1};
                index = index + 2;
            case {'Invert', 'invert'}
                invert = ~invert;
                index = index + 1;
            case {'FlipLR', 'fliplr'}
                flipLR = ~flipLR;
                index = index + 1;
            case 'organizeDepths'
                organizeDepths = ~organizeDepths;
                index = index + 1;
            case {'Verbose', 'verbose'}
                if length(varargin)>index && islogical(varargin{index+1})
                    verbose = varargin{index+1};
                    index = index + 2;
                else
                    verbose = ~verbose;
                    index = index + 1;
                end
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
        Frames(Frames>Config.Frames | Frames<1) = []; % remove requested frames that don't exist
        
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
        
        % Determine indices within file to load
        if Config.Depth>1 % if Config.Depth==1 then frame indices is the same as the file indices
            depthID = idDepth(Config,[],'IndexType',IndexType,'Frames',Frames,'Depths',Depths)';  % determine file indices of frames requested
            numDepths = sum(any(~isnan(depthID),2));
            Frames = sort(depthID(:));                                                          % list of frame indices to load
            Frames(isnan(Frames)) = []; % remove NaN's
        else
            numDepths = 1;
            Frames = Frames';
        end
        numFrames = numel(Frames);

        % Preallocate output
        if info.scanbox_version ~= 1
            Images = zeros(numChannels*HxW*numFrames,1, 'uint16');
        else % version 1
            Images = zeros(numChannels,Config.Width/xavg,Config.Height,1,numFrames, 'uint16');
        end

        % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
        jumps = diff([0;Frames])-1;     % number of frames to jump
        startID = find(jumps);          % indices of jump locations relative to frame indices
        if ~jumps(1)                    % loading very first frame in file
            startID = [1;startID];      % ensures chunk containing first frame is loaded
        end
        jumps = jumps(startID);
        N = diff([startID-1;numFrames]); % number of frames to load in on each loop

        
        % Open File
        info.fid = fopen(SbxFile);
        if(info.fid ~= -1)
            
            % Load Images
            if verbose
                fprintf('Loading   %7.d frame(s) from: %s...', numFrames, SbxFile);
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
                Images = reshape(Images,numChannels,Config.Width,Config.Height,1,numFrames); % reshape vector to usable format
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
            if ~invert % PMTs give a negative signal, so not inverting colormap produces perceptually inverted colormap
                Images = intmax('uint16') - Images;
            end
            
            % Flip across Y axis
            if flipLR
                Images = fliplr(Images);
            end
            
            % Reshape if multiple depths are being loaded
            if organizeDepths && Config.Depth > 1 && numDepths > 1
                Images = depthShape(Images,Frames,depthID','FramesPerDepth',FramesPerDepth);
            end
            
        else
            warning('unable to open file: %s', SbxFile);
            Images = [];
        end
        fclose(info.fid); % close file
        
        
        Config.DimensionOrder = Config.DimensionOrder([3 2 4 1 5]); % update dimensions to that of output as compared to what they are in the file
        Config.size = size(Images); % update size of output
        
        
        if verbose
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