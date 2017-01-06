function Images = depthShape(Images,Frames,depthID)

Dim = size(Images);

if all(diff(Frames)==1) % section is continuous -> simple reshape (assumes acquisition order is 1,2,...,N,1,2,...,N)
    F = find(isnan(depthID'));
    for ii = 1:numel(F) % add in blank frames
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