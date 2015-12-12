function Filename = save2P(Filename,Images,varargin)

Append = true;

Header = ''; % TIF files only

directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Append', 'append'}
                Append = varargin{index+1};
                index = index + 2;
            case {'Header', 'header'}
                Header = varargin{index+1};
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

if ~exist('Images','var') || isempty(Images) % Prompt for file selection
    [Images, directory] = uigetfile({'*.imgs;*.sbx;*.tif'}, 'Select images:', directory, 'MultiSelect','on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
    Images = load2P(Images, 'Type', 'Direct'); % load images
end

% Check saveType
if ~any(strcmp(SaveType,{'tif','sbx'}))
    error('Save type ''%s'' not recognized. Save type has to be either ''tif'' or ''sbx''', SaveType);
end

% Check format
if ~isa(Images(1), 'uint16')
    Images = uint16(Images);
end
numFrames = size(Images, 5);

% Squeeze images
if ndims(Images) == 5
    Images = squeeze(Images(:,:,1,1,:));
end


%% Save images to file

% Determine file type
[~,~,ext] = fileparts(Filename);

% Save images
fprintf('Saving %d frames to file: %s...', numFrames, Filename);
switch ext
    
    case 'tif'
        
        % Open file
        if Append
            tiffObject = Tiff(Filename,'a');
        else
            tiffObject = Tiff(Filename,'w');
        end
        
        % Write frames to file
        for findex=1:numFrames
            if findex~=1
                tiffObject.writeDirectory();
            end
            tiffObject.setTag('ImageLength',size(Images,1));
            tiffObject.setTag('ImageWidth', size(Images,2));
            tiffObject.setTag('Photometric', Tiff.Photometric.MinIsBlack);
            tiffObject.setTag('BitsPerSample', 16);
            tiffObject.setTag('SamplesPerPixel', 1);
            tiffObject.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
            tiffObject.setTag('Software', 'MATLAB');
            if ~isempty(Header)
                tiffObject.setTag('ImageDescription',Header);
            end
            tiffObject.write(Images(:,:,findex));
        end
        
        % Close file
        tiffObject.close();
        
    otherwise % assumes binary file
        
        % Open file
        if Append
            fid = fopen(Filename, 'a');
        else
            fid = fopen(Filename, 'w');
        end
        
        % Permute frames?
        
        % Write frames to file
        for findex=1:numFrames
            fwrite(fid, Images(:,:,findex), 'uint16');
        end
        
        % Close file
        fclose(fid);
        
end %switch saveType

fprintf('\tComplete\n');

