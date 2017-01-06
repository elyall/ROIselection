function outImages = depthShape(Images,Frames,depthID)

Dim = size(Images);

outImages = zeros(Dim(1),Dim(2),size(depthID,2),Dim(4),size(depthID,1),'uint16');

[~,loc] = ismember(Frames,depthID);
[F,D] = ind2sub(size(depthID),loc);
for ii = 1:Dim(5)
    outImages(:,:,D(ii),:,F(ii)) = Images(:,:,:,:,ii);
end