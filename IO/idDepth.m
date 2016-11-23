function Frames = idDepth(Config,varargin)


Frames = inf;
FramesPerDepth = 1;
directory = cd;

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Frames','frames'}
                Frames = varargin{index+1};
                index = index + 2;
            case 'FramesPerDepth'
                FramesPerDepth = varargin{index+1};
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

if ~exist('Config', 'var') || isempty(Config)
    [f,p] = uigetfile({'*.sbx;*.tif'}, 'Choose image file to analyze', directory);
    if isnumeric(f)
        Frames = []; return
    end
    Config = fullfile(p,f);
end


%% Load in config
if ischar(Config)
    Config = load2PConfig(Config);
end


%% Determine frames to load
if ischar(Frames) || (numel(Frames)==1 && Frames == inf)
    Frames = 1:Config.Frames;
elseif Frames(end) == inf
    Frames = [Frames(1:end-2),Frames(end-1):Config.Frames];
end
if Config.Depth>1
    Frames(Frames>ceil(Config.Frames/Config.Depth) | Frames<1) = []; % remove requested frames that don't exist
else
    Frames(Frames>Config.Frames | Frames<1) = []; % remove requested frames that don't exist
end


%% Determine id of frame
PacketSize = FramesPerDepth * Config.Depth;             % number of frames in cycle
F = reshape(1:PacketSize,FramesPerDepth,Config.Depth);	% frame index corresponding to frame 1:FramesPerDepth of each depth
Findex = rem(Frames-1,FramesPerDepth)+1;                % index of each requested frame within FramesPerDepth
Frames = bsxfun(@plus, floor((Frames-1)/FramesPerDepth)'*PacketSize,F(Findex,:)); % frame ID within whole movie; transpose for vectorizing


