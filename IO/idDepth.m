function [FrameIDs,RelativeIndex] = idDepth(totalDepths,totalFrames,varargin)

Depths = [];            % specifies which depths to return absolute frame indices for
Frames = [1,inf];       % specifies which frames the user wants
IndexType = 'absolute'; % 'relative' or 'absolute' -> determines whether index above is the absolute frame indices, or relative frame indices for each depth
FramesPerDepth = 1;     % specifies number of frames taken at given depth before moving on to next depth
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

if ~exist('totalDepths', 'var') || isempty(totalDepths)
    [f,p] = uigetfile({'*.sbx;*.tif'}, 'Choose image file to analyze', directory);
    if isnumeric(f)
        FrameIDs = []; RelativeIndex = []; return
    end
    totalDepths = fullfile(p,f);
end


%% Load in config
if ischar(totalDepths)
    totalDepths = load2PConfig(totalDepths);
end
if isstruct(totalDepths)
    totalFrames = totalDepths.Frames;
    totalDepths = totalDepths.Depth;
end

% In case where there is only one depth, skip the nonsense
if totalDepths == 1
    if Frames(end)==inf
        Frames = [Frames(1:end-2),Frames(end-1):totalFrames]; % load all frames
    end
    FrameIDs = Frames';
    RelativeIndex = FrameIDs;
    return
end


%% Determine id of requested frames

% Determine cycle parameters
CycleSize = FramesPerDepth * totalDepths;                  % number of frames in single cycle
Cycle = reshape(1:CycleSize,FramesPerDepth,totalDepths);   % frame index corresponding to frame 1:FramesPerDepth of each depth

% Determine frame indices
switch IndexType
    case 'relative'

        % Determine relative frames to use
        if rem(totalFrames,CycleSize)                                                                     % last cycle didn't complete
            [RelativeIndicesInLastCycle,~] = find(Cycle<=rem(totalFrames,CycleSize));                     % relative indices of frames in last cycle
        else
            RelativeIndicesInLastCycle = 0;                                                                 % last cycle completed
        end
        maxRelativeIndex = FramesPerDepth*floor(totalFrames/CycleSize)+max(RelativeIndicesInLastCycle);   % maximum relative index in file
        if Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):maxRelativeIndex];                                      % use all frames
        end
        % Frames(Frames>maxRelativeIndex) = [];                                                               % remove requested frames that don't exist
        RelativeIndex = Frames';
       
    case 'absolute'
        
        % Determine absolute frames to use
        if Frames(end) == inf
            Frames = [Frames(1:end-2),Frames(end-1):totalFrames]; % use all frames
        end
        
        % Determine relative indices of frames requested
        [~,ind]=ismember(rem(Frames-1,CycleSize)+1,Cycle);
        [RelativeIndexInCycle,~] = ind2sub(size(Cycle),ind);                            % relative indices of frames in their cycle
        RelativeIndex = FramesPerDepth*floor((Frames-1)/CycleSize)+RelativeIndexInCycle;% relative indices of frames requested
        RelativeIndex = unique(RelativeIndex)';                                          % remove redundant indices
        
end

% Determine frame indices
Findex = rem(RelativeIndex-1,FramesPerDepth)+1;                                                 % index of each requested frame within FramesPerDepth
FrameIDs = bsxfun(@plus, floor((RelativeIndex-1)/FramesPerDepth)*CycleSize,Cycle(Findex,:));   % frame ID within whole movie; transpose for vectorizing

switch IndexType
    case 'relative'
        FrameIDs(FrameIDs>totalFrames) = nan;     % remove frames that don't exist (occurs if last cycle wasn't complete)
    case 'absolute'
        FrameIDs(~ismember(FrameIDs,Frames)) = nan; % remove indices that aren't requested
end


%% Keep only requested depth(s)
if ~isempty(Depths)
    if all(ismember(Depths,1:size(FrameIDs,2)))
        FrameIDs = FrameIDs(:,Depths);
        FrameIDs(all(isnan(FrameIDs),2),:) = []; % remove rows with all frames missing from requested depths
    else
        warning('Depth index requested not contained within depth indices. Returning frame indices for all depths...');
    end
end

