function TifFile = sbx2tif(SbxFile,TifFile,MCdata)

Channel = 1;
Depth = inf;
FrameIndex = [1 inf];

% Memory settings
portionOfMemory = 0.08;     % find 10% or less works best
sizeRAM = 32000000000;      % amount of memory on your computer (UNIX-only)

%% Parse input arguments
if ~exist('SbxFile','var') || isempty(SbxFile)
    [SbxFile,p] = uigetfile({'*.sbx'}, 'Choose sbx file(s) to load', directory, 'MultiSelect', 'on');
    if isnumeric(SbxFile)
        TifFile = [];
        return
    end
    SbxFile = fullfile(p,SbxFile);
elseif ischar(SbxFile)
    [p,f,e] = fileparts(SbxFile);
    if ~isequal(e,'.sbx')
        SbxFile = fullfile(p,[f,'.sbx']);
    end
end

if ~exist('TifFile','var') || isempty(TifFile)
    TifFile = [SbxFile(1:end-3),'tif'];
end

if ~exist('MCdata','var') || isempty(MCdata)
    MCdata = [];
elseif ischar(MCdata)
    MCdata = {MCdata};
end

%% Determine frame size
[Images,loadObj,Config] = load2P(SbxFile,'Frames',2,'Channel',Channel,'Depth',Depth);
sizeFrame = whos('Images');
sizeFrame = sizeFrame.bytes;
if ispc
    mem = memory;
    numFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
else
    numFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
end
numFramesPerLoad = numFramesPerLoad - mod(numFramesPerLoad,Config.Depth); % ensure # of frames loaded in each batch is divisable by the # of depths present (therefore no NaN frames are returned by load2P)
totalFrames = sum([loadObj.files(:).Frames]);


%% Load in motion correction data
if isequal(MCdata,true) % find all align files matching the sbx file
    [p,f,~] = fileparts(SbxFile);
    MCdata = dir(fullfile(p,[f,'*.align']));
    MCdata = fullfile(p,{MCdata(:).name});
end
if iscellstr(MCdata) % load in MCdata from files
    for findex = 1:numel(MCdata)
        temp = load(MCdata{findex},'MCdata','-mat');
        if ~isfield(temp,'MCdata')
            error('File %s does not contain an MCdata variable',MCdata{findex});
        end
        MCdata{findex} = temp.MCdata;
    end
end
if iscell(MCdata) % convert to struct
    MCdata = cat(2,MCdata{:});
end

if ~isempty(MCdata)
    MotionCorrect = true;
else
    MotionCorrect = false;
end

%% Determine frames to process
if FrameIndex(end)==inf
    FrameIndex = cat(2, FrameIndex(1:end-1), FrameIndex(end-1)+1:totalFrames);
end
numFrames = numel(FrameIndex);


%% Cycle through loading, motion correcting, and saving frames in batches
numBatches = max(1,ceil(numFrames/numFramesPerLoad));
if ~iscell(SbxFile)
    fprintf('Saving %d frames from %s to %s in %d batches:\n',numFrames,SbxFile,TifFile,numBatches);
else
    fprintf('Saving %d frames to %s in %d batches:\n',numFrames,TifFile,numBatches);
end
index = 1;
for bindex = 1:numFramesPerLoad:numFrames % direct loading only -> load frames in batches
    fprintf('\t%4.d of %4.d\n',index,numBatches);
    
    % Load current batch
    lastframe = min(bindex+numFramesPerLoad-1,numFrames);
    currentFrames = FrameIndex(bindex:lastframe);
    [Images, loadObj] = load2P(SbxFile,'Frames',currentFrames,'Channel',Channel,'Depth',Depth,'verbose');
    
    % Correct for motion
    if MotionCorrect
        Images = applyMotionCorrection(Images, MCdata, loadObj);
    end
    
    % Re-interleave depths
    [H,W,numZ,numC,N] = size(Images);
    Images = permute(Images,[1,2,4,3,5]);
    Images = reshape(Images,H,W,numC,numZ*N);
    if lastframe==numFrames % remove all NaN frames
        N = numZ*N-numel(currentFrames);
        Images(:,:,:,end-N+1:end) = [];
    end
    
    % Save to file
    if bindex==1
        save2P(TifFile,Images,'verbose');
    else
        save2P(TifFile,Images,'append','verbose');
    end
    
    index = index + 1;
end

fprintf('Complete\n');
