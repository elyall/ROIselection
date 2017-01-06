function FrameIDs = idDepth_simple(Config,varargin)

Depths = [];
Frames = [1,inf];
IndexType = 'absolute'; % 'relative' or 'absolute'
directory = cd;

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Depth','Depths'}
                Depths = varargin{index+1};
                index = index + 2;
            case {'Frame','Frames'}
                Frames = varargin{index+1};
                index = index + 2;
            case 'IndexType'
                IndexType = varargin{index+1};
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

% In case where there is only one depth, skip the nonsense
if Config.Depth == 1
    if Frames(end)==inf
        Frames = [Frames(1:end-2),Frames(end-1):Config.Frames]; % load all frames
    end
    Frames = Frames';
    return
end


%% Determine id of requested frames


% Determine frame indices
switch IndexType
    case 'relative'

        % Determine relative frames to use
        maxRelativeIndex = floor((Config.Frames-1)/Config.Depth)+1; % maximum relative index in file
        if Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):maxRelativeIndex];  % use all frames
        end
        FrameIDs = bsxfun(@plus,(Frames-1)'*Config.Depth,1:Config.Depth);
        FrameIDs(FrameIDs>Config.Frames) = nan; % remove frames that don't exist (occurs if last cycle wasn't complete)
       
    case 'absolute'
        
        % Determine absolute frames to use
        if Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):Config.Frames]; % use all frames
        end
        
        % Determine relative indices of frames requested
        RelativeIndex = floor((Frames-1)/Config.Depth)+1;
        RelativeIndex = unique(RelativeIndex)';
        FrameIDs = bsxfun(@plus,(RelativeIndex-1)*Config.Depth,1:Config.Depth);
        FrameIDs(~ismember(FrameIDs,Frames)) = nan; % remove indices that aren't requested

end


%% Keep only requested depth(s)
if ~isempty(Depths)
    if all(ismember(Depths,1:size(FrameIDs,2)))
        FrameIDs = FrameIDs(:,Depths);
    else
        warning('Depth index requested not contained within depth indices. Returning frame indices for all depths...');
    end
end

