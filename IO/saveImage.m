function Filename = saveImage(Filename,Images,varargin)

Append = false;
Scale = true;

directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Append', 'append'}
                Append = varargin{index+1};
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

if ~exist('Filename','var') || isempty(Filename) % No filename input
    [Filename, directory] = uiputfile({'*.tif;*.sbx'}, 'Save images as:', directory);
    if isnumeric(Filename)
        return
    end
    Filename = fullfile(directory,Filename);
end
[~,~,ext] = fileparts(Filename); % Determine file type

if ~exist('Images','var') || isempty(Images) % Prompt for file selection
    [Images, directory] = uigetfile({'*.tif'}, 'Select image file:', directory);
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
end


%% Load and format images

% Load images
if ischar(Images)
    info = imfinfo(Images,'tif');
    tif=Tiff(Images,'r');
    Images = zeros(info(1).Height, info(1).Width, 1, numel(info));
    for findex = 1:numel(info)
        if findex ~= 1
            tif.nextDirectory;
        end
        Images(:,:,1,findex)=tif.read();
    end
    tif.close;
end

% Format images
if ndims(Images)==3 && size(Images,3)~=3
    Images = permute(Images, [1,2,4,3]);
end
[H,W,~,numFrames] = size(Images);

% Scale images
if Scale
    if isa(Images,'double')
        MaxInt = double(intmax('uint32'));
    elseif isa(Images,'single')
        MaxInt = single(intmax('uint16'));
    else
        MaxInt = intmax(class(Images));
    end
    Images = Images - min(Images(:)) + 1;    % shift bottom to 1
    Images = Images*MaxInt/max(Images(:));   % scale top to maximum value of bit depth
end

% Determine bit depth
switch class(Images)
    case {'double','uint32','int32'}
        bitDepth = 32;
        if isa(Images, 'double') && strcmp(ext,'.tif')
            Images = uint32(Images);
        end 
    case {'single','uint16','int16'}
        bitDepth = 16;
    case {'uint8','int8'}
        bitDepth = 8;
end


%% Save images to file

% Save images
if Append
    fprintf('Appending\t%d\tframe(s) to file:\t%s...', numFrames, Filename);
else
    fprintf('Writing\t%d\tframe(s) to file:\t%s...', numFrames, Filename);
end

switch ext
    
    case '.tif'
        
        % Open file
        if Append
            tiffObject = Tiff(Filename,'a');
        else
            tiffObject = Tiff(Filename,'w');
        end
        
        % Write frames to file
        for findex = 1:numFrames
            if findex~=1
                tiffObject.writeDirectory();
            end
            tiffObject.setTag('ImageLength',H);
            tiffObject.setTag('ImageWidth', W);
            tiffObject.setTag('Photometric', Tiff.Photometric.MinIsBlack);
            tiffObject.setTag('BitsPerSample', bitDepth);
            tiffObject.setTag('SamplesPerPixel', 1);
            tiffObject.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
            tiffObject.setTag('Software', 'MATLAB');
            tiffObject.write(Images(:,:,:,findex));
        end
        
        % Close file
        tiffObject.close();
        
end %switch ext

fprintf('\tComplete\n');

