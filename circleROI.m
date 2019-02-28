function [ROImask, Position] = circleROI(Image,Position)


hF = figure;

hA = axes('Units','Normalized','Position',[0,0,1,.9]);
imagesc(Image); hold on; axis off;

hB = [];
hB(1) = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Add ROI',...
    'Parent',               hF,...
    'Units',                'normalized',...
    'Position',             [0,.9,.33,.1],...
    'Callback',             @(hObject,eventdata)AddROI(hObject,eventdata,guidata(hObject)));
hB(2) = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Save ROI',...
    'Parent',               hF,...
    'Units',                'normalized',...
    'Position',             [.33,.9,.34,.1],...
    'Enable',               'off',...
    'Callback',             @(hObject,eventdata)SaveROI(hObject,eventdata,guidata(hObject)));
hB(3) = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Delete ROI',...
    'Parent',               hF,...
    'Units',                'normalized',...
    'Position',             [.67,.9,.33,.1],...
    'Enable',               'off',...
    'Callback',             @(hObject,eventdata)DeleteROI(hObject,eventdata,guidata(hObject)));

h = [];
dict = [];
ROImask = [];
if exist('Position','var') && ~isempty(Position)
    PlotInputROIs(Position)
else
    Position = [];
end
fcnE = makeConstrainToRectFcn('imellipse',get(gca,'XLim'),get(gca,'YLim'));
% fcnP = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));

waitfor(hF);

    function PlotInputROIs(Position)
        for r = 1:size(Position,1)
            h = imellipse(hA,Position(r,:));
            ROImask = cat(3,ROImask,createMask(h));
            vertices = getVertices(h);
            delete(h);
            dict = cat(1,dict,datetime);
            patch(vertices(:,1),vertices(:,2),[1,0,0],...
                'Parent',    hA,...
                'FaceAlpha', 0,...
                'EdgeColor', [1,0,0],...
                'LineWidth', 1,...
                'UserData',  dict(end),...
                'ButtonDownFcn',@(hObject,eventdata)EditROI(hObject,eventdata,guidata(hObject)));
        end
    end

    function AddROI(hObject,eventdata,gd)
        h = imellipse(hA);
        setPositionConstraintFcn(h, fcnE);
        % h = impoly(gca,'Closed',1);
        % setPositionConstraintFcn(h, fcnP);
        set(hB(1),'Enable','off');
        set(hB(2:3),'Enable','on');
    end

    function SaveROI(hObject,eventdata,gd)
        ROImask = cat(3,ROImask,createMask(h)); % extract mask
        Position = cat(1,Position,getPosition(h));
        vertices = getVertices(h);
        delete(h);
        dict = cat(2,dict,datetime);
        patch(vertices(:,1),vertices(:,2),[1,0,0],...
            'Parent',    hA,...
            'FaceAlpha', 0,...
            'EdgeColor', [1,0,0],...
            'LineWidth', 1,...
            'UserData',  dict(end),...
            'ButtonDownFcn',@(hObject,eventdata)EditROI(hObject,eventdata,guidata(hObject)));
        set(hB(1),'Enable','on');
        set(hB(2:3),'Enable','off');
    end

    function EditROI(hObject,eventdata,gd)
        id = find(dict==hObject.UserData);
        h = imellipse(hA,Position(id,:));
        ROImask(:,:,id) = [];
        Position(id,:) = [];
        dict(id) = [];
        delete(hObject);
        setPositionConstraintFcn(h, fcnE);
        set(hB(1),'Enable','off');
        set(hB(2:3),'Enable','on');
    end

    function DeleteROI(hObject,eventdata,gd)
        delete(h);
        set(hB(1),'Enable','on');
        set(hB(2:3),'Enable','off');
    end

end