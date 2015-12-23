function Filename = save2P(Filename,Images,varargin)

Append = false;
Header = ''; % tif: header, bin: mat-file
Invert_binary = true; % bin-only: invert colormap before saving

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
            case 'Invert'
                Invert_binary = varargin{index+1};
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

% Check format
if ~isa(Images(1), 'uint16')
    Images = uint16(Images);
end
[~,~,numZ,numC,numFrames] = size(Images);


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
    
    case 'tif'
        
        % Open file
        if Append
            tiffObject = Tiff(Filename,'a');
        else
            tiffObject = Tiff(Filename,'w');
        end
        
        % Write frames to file
        for dindex = 1:numZ
            for cindex = 1:numC
                for findex = 1:numFrames
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
                    tiffObject.write(Images(:,:,dindex,cindex,findex));
                end
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
        Images = permute(Images, [4, 2, 1, 3, 5]); % change to [channel, width, height, depth, frame]
        
        % Invert colormap
        if Invert_binary
            Images = intmax('uint16') - Images;
        end
        
        % Write frames to file
        fwrite(fid, Images, 'uint16');
        
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

