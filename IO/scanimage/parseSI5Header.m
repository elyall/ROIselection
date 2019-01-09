function header = parseSI5Header(input)

sIndex=strfind(input,'SI.');
vs =5.2;

if isempty(sIndex)
    sIndex=strfind(input,'scanimage.');
    vs = 5;
end

if vs==5
    for i=1:size(sIndex,2)-1
        
        [field, remainder ] = strtok(input(sIndex(i):sIndex(i+1)-1));
        [tk, val] = strtok(remainder);
        field= strrep(field,'.','_'); %remove . from names if present
        val = val(2:end-1);
        try
            eval(['header.scanimage.(field(15:end))=' val ';']);
        catch %catch weird variables and translate to strings <>
            header.scanimage.(field(15:end))=val;
        end
    end
    [field, remainder ] = strtok(input(sIndex(i+1):end));
    [tk, val] = strtok(remainder);
    field= strrep(field,'.','_'); %remove . from names if present
    val2 = strtok(val);
    try
        eval(['header.scanimage.(field(15:end))=' val2 ';']);
    catch
        header.scanimage.(field(15:end))=val;
    end
    
elseif vs ==5.2 %should be largely similar to vs 5 just a few changes
    for i=1:size(sIndex,2)-1
        
        [field, remainder ] = strtok(input(sIndex(i):sIndex(i+1)-1));
        [tk, val] = strtok(remainder);
        %field= strrep(field,'.','_'); %remove . from names if present
        val = val(2:end-1);
        try
            eval(['header.SI.' field(4:end) ' =' val ';']);
        catch %catch weird variables and translate to strings <>
            header.SI.(field(4:end))=val;
        end
    end
    [field, remainder ] = strtok(input(sIndex(i+1):end));
    [tk, val] = strtok(remainder);
    %field= strrep(field,'.','_'); %remove . from names if present
    val2 = strtok(val);
    try
        eval(['header.SI.' field(4:end) ' = ' val2 ';']);
    catch
        header.SI.(field(4:end))=val;
    end
end