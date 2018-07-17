function Images = depthShape(Images,Frames,depthID,varargin)
%DEPTHSHAPE De-interleaves frames of different depths
%   IMAGES = depthShape(IMAGES,FRAMES,DEPTHID) reshapes IMAGES of dimension
%   [H x W x 1 x C x F*D] into [H x W x D x C x F]. FRAMES is a vector
%   specifying the absolute frame index of each frame in IMAGES. DEPTHID is
%   a matrix returned by idDepth for absolute indices FRAMES specifying the
%   absolute frame index of each relative frame for each depth (see
%   idDepth).
%
%   IMAGES = depthShape(...,'FramesPerDepth',FPD) specifies how many frames
%   were taken at each depth before moving to the next depth.
%


% Default parameters that can be adjusted
FramesPerDepth = 1; % specifies number of frames taken at given depth before moving on to next depth

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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

% For sampling multiple frames at each depth -> reorder so that depths are
% interleaved every frame (simplifies reordering)
if numel(FramesPerDepth)>1 && numel(unique(FramesPerDepth))==1 && FramesPerDepth(1)>1
    order = depthID';
    order = order(:);
    Images = Images(:,:,:,:,order);
    FramesPerDepth = 1;
end

%% De-interleave different depths
Dim = size(Images);
if all(FramesPerDepth==1) % section is continuous -> simple reshape 
    F = find(isnan(depthID'));  % locate nonexistent frames
    for ii = 1:numel(F)         % add in blank frames
        if F(ii) == 1
            Images = cat(5,zeros(Dim(1:4),'uint16'),Images);
        elseif F(ii) > size(Images,5)
            Images = cat(5,Images,zeros(Dim(1:4),'uint16'));
        else
            Images = cat(5,Images(:,:,:,:,1:F(ii)-1),zeros(Dim(1:4),'uint16'),Images(:,:,:,:,F(ii):end));
        end
    end
    Images = reshape(Images, [Dim(1:2),size(depthID,2),Dim(4),size(depthID,1)]); % reshape data
    
else % data is discontinuous
    out = zeros(Dim(1),Dim(2),size(depthID,2),Dim(4),size(depthID,1),'uint16');	% initialize output
    [~,loc] = ismember(Frames,depthID);                                         % determine order
    [F,D] = ind2sub(size(depthID),loc);                                         % change order to frame and depth indices
    for ii = 1:Dim(5)
        out(:,:,D(ii),:,F(ii)) = Images(:,:,:,:,ii);                            % assign frames to proper location
    end
    Images = out;
end

