function ROImask = circleROI(Image)


hF = figure;

hA = axes('Units','Normalized','Position',[0,0,1,.9]);
imagesc(Image);

hB = uicontrol(...
    'Style',                'togglebutton',...
    'String',               'Add ROI',...
    'Parent',               hF,...
    'Units',                'normalized',...
    'Position',             [0,.9,1,.1],...
    'Callback',             @(hObject,eventdata)ROI(hObject,eventdata,guidata(hObject)));

h = [];
ROImask = [];

waitfor(hF);

    function ROI(hObject,eventdata,gd)
        if hObject.Value
            fcnE = makeConstrainToRectFcn('imellipse',get(gca,'XLim'),get(gca,'YLim'));
            h = imellipse(gca);
            setPositionConstraintFcn(h, fcnE);
            % fcnP = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
            % h = impoly(gca,'Closed',1);
            % setPositionConstraintFcn(h, fcnP);
            
            hObject.String = 'Save ROI';
        else
            ROImask = cat(3,ROImask,createMask(h)); % extract mask
            delete(h);
            hObject.String = 'Add ROI';
        end
    end

end