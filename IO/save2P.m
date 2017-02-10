function Filename = save2P(Filename,Images,varargin)
%SAVE2P Saves images to .tif or .sbx file
%   save2P() prompts user to select a filename to save to and select a .sbx
%   or .tif file to load data from.
%
%   save2P(FILENAME,IMAGES) determines what the name of the file being
%   saved to and the images being saved to it.
%
%   save2P(...,'Header',HEADER) saves the metadata HEADER. For tif files
%   HEADER should be a string, for sbx files HEADER should be a struct.
%   HEADER should be properly formatted if user wants to load the file
%   without issue. Leave empty if don't want metadata saved (default='')
%
%   save2P(...,'append') appends the IMAGES to the end of the file rather
%   than overwriting the file.
%
%   save2P(...,'class',CLASS) string that sets IMAGES to be of class CLASS.
%   (default='').
%
%   save2P(...,'invert') inverts the colormap of the data before saving it
%   to file.
%

% Default parameters that can be adjusted
Append = false; % booleon determing whether data is written to a new file (potentially overwriting existing file) or appended to existing one
Header = '';    % metadata -> tif: string saved as header for each frame, bin: struct saved to corresponding mat-file
invert = true;  % booleon determing whether to invert the colormap before saving
Class = '';     % string specifying class to save the images as

% Placeholders
directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Append', 'append'}
                Append = ~Append;
                index = index + 1;
            case {'Header', 'header'}
                Header = varargin{index+1};
                index = index + 2;
            case {'Invert','invert'}
                invert = ~invert;
                index = index + 1;
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
end


%% Load images
if ischar(Images) || iscellstr(Images)
    Images = load2P(Images, 'Type', 'Direct'); % load images
end

% Set format
if ~isempty(Class)
    if ~isa(Images, Class)
        switch Class
            case 'uint8'
                Images = uint8(Images);
            case 'uint16'
                Images = uint16(Images);
            case 'uint32'
                Images = uint32(Images);
            case 'uint64'
                Images = uint64(Images);
            case 'single'
                Images = single(Images);
            case 'double'
                Images = double(Images);
        end
    end
end
Class = class(Images);

% Invert colormap
if invert
    Images = intmax(Class) - Images;
end

% Re-interleave depths
[H,W,numZ,numC,numFrames] = size(Images);
Images = permute(Images,[1,2,4,3,5]);
Images = reshape(Images,H,W,numC,numZ*numFrames);


%% Save images to file

% Determine file type
[~,~,ext] = fileparts(Filename);

% Save images
if Append
    fprintf('Appending\t%d\tframes to file:\t%s...', numFrames, Filename);
else
    fprintf('Writing\t%d\tframes to file:\t%s...', numFrames, Filename);
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
        for findex = 1:numFrames*numZ
            for cindex = 1:numC
                if findex~=1 || cindex~=1
                    tiffObject.writeDirectory();
                end
                tiffObject.setTag('ImageLength',H);
                tiffObject.setTag('ImageWidth',W);
                tiffObject.setTag('Photometric', Tiff.Photometric.MinIsBlack);
                tiffObject.setTag('BitsPerSample', 16);
                tiffObject.setTag('SamplesPerPixel', 1);
                tiffObject.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
                tiffObject.setTag('Software', 'MATLAB');
                if ~isempty(Header)
                    tiffObject.setTag('ImageDescription',Header);
                end
                tiffObject.write(Images(:,:,cindex,findex));
            end
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
        
        % Permute frames
        Images = permute(Images, [3, 2, 1, 4]); % change to [channel, width, height, frame]
        
        % Write frames to file
        fwrite(fid, Images, Class);
        
        % Close file
        fclose(fid);
        
        % Save header file
        if ~isempty(Header)
            [p,f,~] = fileparts(Filename);
            HeaderFile = fullfile(p,[f '.mat']);
            info = Header;
            save(HeaderFile, 'info', '-mat', '-v7.3');
        end
        
end %switch ext

fprintf('\tComplete\n');

