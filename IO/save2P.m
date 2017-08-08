function Filename = save2P(Filename,Images,varargin)
%SAVE2P Saves images to .tif, .avi, or binary file
%   FILENAME = save2P() prompts user to select a FILENAME to save to and
%   select a .sbx or .tif file to load data from.
%
%   save2P(FILENAME,IMAGES) sets the name of the file to save to and the
%   images to save to it. Images can be 2D [Y,X], 3D [Y,X,F], 4D [Y,X,C,F],
%   or 5D [Y,X,Z,C,F]. FILENAME can specify a '.tif' file, a '.avi' file,
%   or otherwise saves the data as a binary file.
%
%   save2P(...,'Header',HEADER) saves the metadata HEADER. For tif files
%   HEADER should be a string, for sbx files HEADER should be a struct.
%   HEADER should be properly formatted if user wants to load the file
%   without issue. Leave empty if don't want metadata saved. (default='')
%
%   save2P(...,'append') appends the IMAGES to the end of the file rather
%   than overwriting the file.
%
%   save2P(...,'CLim',CLIM) either a 1x2 vector specifying the color limits
%   to apply to IMAGES, or True to scale each frame to use the full
%   colormap. (default=[])
%
%   save2P(...,'class',CLASS) string that sets IMAGES to be of class CLASS.
%   (default='')
%
%   save2P(...,'invert') inverts the colormap of the data before saving it
%   to file.
%
%   save2P(...,'frameRate',FRAMERATE) specifies the frame rate of the AVI
%   file. (default=15)
%

% Default parameters that can be adjusted
Header = '';     % metadata -> tif: string saved as header for each frame, bin: struct saved to corresponding mat-file
Append = false;  % booleon determing whether data is written to a new file (potentially overwriting existing file) or appended to existing one
CLim = [];       % vector of length 2 specifying the limits of the colormap
MaxOutValue = [];% scalar specifying what to set the output max to be
Class = '';      % string specifying class to save the images as
invert = false;  % booleon determing whether to invert the colormap before saving
frameRate = 15;  % avi only: scalar specifying frameRate of video

% Placeholders
directory = cd; % default directory when prompting user to select a file

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
            case {'Class','class'}
                Class = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
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
end


%% Load & adjust images
if ischar(Images) || iscellstr(Images)
    if all(contains(Images,{'.tif','.tiff','.sbx'}))
        Images = load2P(Images, 'Type', 'Direct'); % load images
    elseif all(contains(Images,'.align'))
        Files = Images;
        Images = [];
        for findex = 1:numel(Files)
            load(Files{findex},'m','-mat')
            Images = cat(3,Images,m);
        end
    end
end

% Format Images
if ndims(Images)==3     % add channel dimension
    Images = permute(Images,[1,2,4,3]);
elseif ndims(Images)==5 % re-interleave depths
    [H,W,numZ,numC,numFrames] = size(Images);
    Images = permute(Images,[1,2,4,3,5]);
    Images = reshape(Images,H,W,numC,numZ*numFrames);
elseif ndims(Images)>5
    error('Images must have less than 5 dimensions!');
end
[H,W,numC,numFrames] = size(Images);

% Set color limits
if isequal(CLim,true)
%     CLim = [min(Images(:)),max(Images(:))];
    CLim = prctile(Images(:), [.1,99.9]);
end
if ~isempty(CLim)
    if isempty(Class)
        Class = class(Images);
    end
    Images = double(Images);
    if isequal(CLim,true) % imagesc each frame
        for findex = 1:numFrames
            for cindex = 1:numC
                temp = Images(:,:,cindex,findex);
                Images(:,:,cindex,findex) = (temp-min(temp(:)))./range(temp(:));
            end
        end
    else
        Images = (Images - CLim(1))./diff(CLim);
    end
    Images(Images>1) = 1; % with max at 1, rectify all larger values
    Images(Images<0) = 0; % with min at 0, rectify all smaller values
    if isempty(MaxOutValue)
        switch Class
            case {'single','double'}
                Images = (double(realmax(Class))-1)*Images+1;
            otherwise
                Images = (double(intmax(Class))-1)*Images+1;
        end
    else
        Images = (double(MaxOutValue)-1)*Images+1;
    end
end

% Set format
if ~isempty(Class) && ~isa(Images, Class)
    Images = cast(Images, Class); % cast to desired format
end
Class = class(Images);

% Invert colormap
if invert
    switch Class
        case {'single','double'}
            Images = realmax(Class) - Images;
        otherwise
            Images = intmax(Class) - Images;
    end
end


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
        
        % Ensure data works with TiffLib
        if isa(Images,'single')
            Images = uint32(Images);
            Class = 'uint32';
        elseif isa(Images,'uint64') || isa(Images,'int64')
            Images = double(Images);
            Class = 'double';
        end

        % Create Tiff metadata
        meta.Photometric = Tiff.Photometric.MinIsBlack;
        meta.ImageLength = H;
        meta.ImageWidth = W;
        meta.RowsPerStrip = H;
        switch Class
            case 'logical'
                meta.BitsPerSample = 1;
            case {'int8','uint8'}
                meta.BitsPerSample = 8;
            case {'int16','uint16'}
                meta.BitsPerSample = 16;
            case {'int32','uint32'}
                meta.BitsPerSample = 32;
            case 'double'
                meta.BitsPerSample = 64;
        end
        meta.Compression = Tiff.Compression.None;
        switch Class
            case {'logical','uint8','uint16','uint32'}
                meta.SampleFormat = Tiff.SampleFormat.UInt;
            case {'int8','int16','int32'}
                meta.SampleFormat = Tiff.SampleFormat.Int;
            case 'double'
                meta.SampleFormat = Tiff.SampleFormat.IEEEFP;
        end
        meta.SamplesPerPixel = 1;
        meta.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        meta.Software = 'MATLAB';
        if ~isempty(Header)
            meta.ImageDescription = Header;
        end
        
        % Open file
        if Append
            tiffObject = Tiff(Filename,'a');
        else
            tiffObject = Tiff(Filename,'w');
        end
                
        % Write frames to file
        for findex = 1:numFrames
            for cindex = 1:numC
                if findex~=1 || cindex~=1
                    tiffObject.writeDirectory();
                end
                tiffObject.setTag(meta);                     % set metadata for current directory
                tiffObject.write(Images(:,:,cindex,findex)); % save current image
            end
        end
        
        % Close file
        tiffObject.close();
        
    case '.avi'
        
        % Ensure data works with VideoWriter
        if isa(Images,'uint16')
            Images = single(Images);
            Class = 'single';
        elseif isa(Images,'uint32')
            Images = single(Images);
            Class = 'single';
        elseif isa(Images,'uint64')
            Images = double(Images);
            Class = 'double';
        end
        
        % Set to be between 0 and 1
        if any(strcmp(Class,{'single','double'})) && (max(Images(:))>1 || min(Images(:))<0)
            Images = Images - min(Images(:)); % set min to 0
            Images = Images/max(Images(:));   % set max to 1
        end
        
        % Write video
        vidObj = VideoWriter(Filename,'Motion JPEG AVI');   % initiate file
        set(vidObj, 'FrameRate', frameRate);                % set frame rate
        open(vidObj);                                       % open file
        for findex = 1:numFrames
            for cindex = 1:numC
                writeVideo(vidObj, Images(:,:,cindex,findex)); % write frame
            end
        end
        close(vidObj);                                      % close file
    
%     case 'hdf5'
%         
%         sizY = [Y,X,F];
%         bin_width = round(512*512/prod(sizY)*4e3);
%         nd = length(sizY)-1;
%         sizY = sizY(1:nd);
%         h5create(output_filename,'/mov',[sizY,Inf],'Chunksize',[sizY,bin_width],'Datatype',data_type);
%         for t = 1:bin_width:T
%         h5write(output_filename,'/mov',Ytm,[ones(1,nd),t],[sizY(1:nd),size(Ytm,3)]);
%         end
        
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
        
        % Save header file (sbx format)
        if ~isempty(Header)
            [p,f,~] = fileparts(Filename);
            HeaderFile = fullfile(p,[f '.mat']);
            info = Header;
            save(HeaderFile, 'info', '-mat');
        end
        
end %switch ext

fprintf('\tComplete\n');

