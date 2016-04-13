function varargout = intelligentScissor(varargin)
% INTELLIGENTSCISSOR MATLAB code for intelligentScissor.fig
%      INTELLIGENTSCISSOR, by itself, creates a new INTELLIGENTSCISSOR or raises the existing
%      singleton*.
%
%      H = INTELLIGENTSCISSOR returns the handle to a new INTELLIGENTSCISSOR or the handle to
%      the existing singleton*.
%
%      INTELLIGENTSCISSOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTELLIGENTSCISSOR.M with the given input arguments.
%
%      INTELLIGENTSCISSOR('Property','Value',...) creates a new INTELLIGENTSCISSOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before intelligentScissor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to intelligentScissor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help intelligentScissor

% Last Modified by GUIDE v2.5 03-Mar-2016 00:24:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @intelligentScissor_OpeningFcn, ...
                   'gui_OutputFcn',  @intelligentScissor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before intelligentScissor is made visible.
function intelligentScissor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to intelligentScissor (see VARARGIN)

% Choose default command line output for intelligentScissor
handles.output = hObject;
handles.scale = 1;
handles.mode = 'original';
set(handles.imageView,'xtick',[]);
set(handles.imageView,'ytick',[]);
% handles.node = struct('pixel',{},'linkCost',zeros(1,8),...
%     'state','INITIAL','totalCost',0,'PrevNode',zeros(1,2),...
%     'col',0,'row',0);
%set(handles.imageView,'Position', [25, 100, 375, 400]);
addpath('export_fig');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes intelligentScissor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = intelligentScissor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in startTracing.
function startTracing_Callback(hObject, eventdata, handles)
% hObject    handle to startTracing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curFig = get(0, 'CurrentFigure');
figure(curFig);

[handles.height, handles.width, handles.dim] = size(handles.current);
if handles.dim == 3
    Img = rgb2gray(handles.current);
else
    Img = handles.current;
end
Img = imgaussfilt(Img, 5);
try
    set(curFig, 'WindowButtonDownFcn', @liveWire, 'WindowButtonMotionFcn', @mouseMove2, ...
        'KeyPressFcn', @keyFuncs1,'DoubleBuffer', 'on'); 
catch
    set(curFig, 'WindowButtonDownFcn', @liveWire, 'WindowButtonMotionFcn', @mouseMove2, ...
        'KeyPressFcn', @keyFuncs1, 'DoubleBuffer', 'on');
end

handles.contours(handles.maskInd,1) = line('Parent', gca, 'XData', [], 'YData', [], 'Clipping', 'off', ...
    'Color', 'g', 'LineStyle', '-', 'LineWidth', 1.5);
handles.contours(handles.maskInd,2) = line('Parent', gca, 'XData', [], 'YData', [], 'Clipping', 'off', ...
    'Color', 'y', 'LineStyle', ':', 'LineWidth', 1.5);

handles.L = getCost(Img); %local costs
% path map with lowest cost
handles.xPathMap = zeros(size(handles.current), 'int8');
handles.yPathMap = zeros(size(handles.current), 'int8'); 
% path coordinates
handles.XpathPos = []; 
handles.YpathPos = []; 
%seed point list
handles.seedList = zeros(200, 1);
handles.seedIndx  = 0;
guidata(hObject, handles);


function viewMode_Callback(hObject, eventdata, handles)
% hObject    handle to viewMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns viewMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from viewMode
% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
switch str{val}
    case 'Original Image' 
        handles.mode = 'original';
    case 'Gradient Map'
        handles.mode = 'gradient';    
end
viewModeSet(hObject, eventdata, handles);
% Save the handles structure.
guidata(hObject,handles)

function viewMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [filename1,filepath1]=uigetfile({'*.*','All Files'},'Select File');
  if filename1
      cd(filepath1);
      handles.g = imread(filename1);
      handles.current = handles.g;
      axes(handles.imageView);
      handles.scale = 1;
      [handles.height, handles.width, handles.dim] = size(handles.current);
      set(handles.imageSize,'String',['Image Size: ',int2str(handles.width),'*',int2str(handles.height),...
          char(10),'Image Dimension: ',int2str(handles.dim),char(10),'Scale: ',num2str(handles.scale)]);
      hold all;
      xlim([-1 2.5*handles.width]);
      ylim([-1 2.5*handles.height]);
      delete(get(gca,'Children'));
      imshow(handles.g);
      handles.masks = [];
      handles.maskInd = 1;
      handles.selectedMaskInd = 0;
      handles.coordinates = cell(0);
      handles.contours = [];
      handles.C1 = [];
      guidata(hObject, handles);
  end

function saveContour_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file_name,filepath] = uiputfile({'*.png;*.tif;*.jpg;*.gif','All Image Files';...
          '*.*','All Files' },'Save Image','contour.png');
if file_name
    cd(filepath);
    h = findobj(gcf,'type','axes');
    f1 = figure('visible','off'); 
    copyobj(h,f1); 
    f1a = findobj(f1,'type','axes');
    axis image 
    export_fig(fullfile(file_name),'-native');
%     %screenshot the axes frame
%     F = getframe(f1a);
%     Image = frame2im(F);    
%     imwrite(Image, filename);
    clf(f1);
end

function saveMask_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.selectedMaskInd == 0
    h = msgbox('Select a mask first!');
else
    f1 = figure('visible','on'); 
    M = handles.masks(:,:,handles.selectedMaskInd);
    imshow(M);
    [file_name,filepath] = uiputfile({'*.mat','All Files'},'Save Mask','mask.mat');
    if file_name
        cd(filepath);
        save(file_name,'M');
        axis image 
    end
end

% --------------------------------------------------------------------
function saveCuttedImg_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveCuttedImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.selectedMaskInd == 0
    h = msgbox('Select a mask first!');
else
    [file_name,filepath] = uiputfile({'*.png;*.tif;*.jpg;*.gif','All Image Files';...
          '*.*','All Files' },'Save Cutted Image','cuttedImg.png');
    if file_name
     cd(filepath);
     M = handles.masks(:,:,handles.selectedMaskInd);
     Inew = bsxfun(@times, handles.current, cast(M, class(handles.current)));
     image(Inew,'alphadata',M);
     imwrite(Inew,file_name,'png','Transparency',[0 0 0]);
    end
end


% function mouseMove (hObject, eventdata, handles)
% handles = guidata(hObject);
% handles.C = get (gca, 'CurrentPoint');
% set(handles.pos,'String',['Cursor Position', char(10),'X: ',num2str(handles.C(1,1)),...
%     char(10),'Y: ',num2str(handles.C(1,2))]);
% guidata(hObject, handles);




function zoomOut_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to zoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
if handles.scale < 2
    handles.scale = handles.scale * 1.25;
    handles.current = imresize(handles.g,handles.scale);
    delete(get(gca,'Children'));
    imshow(handles.current);
    [handles.height, handles.width, handles.dim] = size(handles.current);
    set(handles.imageSize,'String',['Image Size: ',int2str(handles.width),'*',int2str(handles.height),...
      char(10),'Image Dimension: ',int2str(handles.dim),char(10),'Scale: ',num2str(handles.scale)]);
    viewModeSet(hObject, eventdata, handles);
    guidata(hObject, handles);
end

function zoomIn_ClickedCallback(hObject, eventdata, handles)
if handles.scale > 0.32
    handles.scale = handles.scale * 0.8;
    handles.current = imresize(handles.g,handles.scale);
    delete(get(gca,'Children'));
    imshow(handles.current);
    [handles.height, handles.width, handles.dim] = size(handles.current);
    set(handles.imageSize,'String',['Image Size: ',int2str(handles.width),'*',int2str(handles.height),...
      char(10),'Image Dimension: ',int2str(handles.dim),char(10),'Scale: ',num2str(handles.scale)]);
    viewModeSet(hObject, eventdata, handles);
    guidata(hObject, handles);
end

function viewModeSet(hObject, eventdata, handles)
%handles = guidata(hObject);
switch handles.mode;
case 'original' 
    handles.current = imresize(handles.g,handles.scale);
    delete(get(gca,'Children'));
    imshow(handles.current);
case 'gradient' 
    img = rgb2gray(imresize(handles.g,handles.scale));
    prewittKernelV = [-1 -1 -1; 0 0 0; 1 1 1];
    outputImageV = imfilter(img, prewittKernelV);
    prewittKernelH = [-1 0 1; -1 0 1; -1 0 1];
    outputImageH = imfilter(img, prewittKernelH);
    handles.current = uint8(sqrt(single(outputImageH).^2 + single(outputImageV).^2));
    delete(get(gca,'Children'));
    imshow(handles.current);
end
guidata(hObject, handles);


function [goodX, goodY] = seedSnapping(hObject, eventdata, handles, Xpos, Ypos, approxR)
handles = guidata(hObject);
if Xpos < 0.5 ||Ypos < 0.5||Xpos > handles.width||Ypos > handles.height, return, end
Xpos = int16(Xpos);
Ypos = int16(Ypos);
xLowerLim = max(1, Xpos - approxR);
xUpperLim = min(size(handles.L, 2), Xpos + approxR);
yLowerLim = max(1, Ypos - approxR);
yUpperLim = min(size(handles.L, 1), Ypos + approxR);
minIndX = Xpos;
minIndY = Ypos;
for i = xLowerLim: xUpperLim
    for j = yLowerLim: yUpperLim
        if handles.L(i,j) < handles.L(minIndX, minIndY);
            minIndX = i;
            minIndY = j;
        end
    end
end
goodX = minIndX;
goodY = minIndY;
        
        
function liveWire(hObject, eventdata, handles)
global tCost seedX seedY;

handles = guidata(hObject);
C = get (gca, 'CurrentPoint');
Xpos = C(1,1);
Ypos = C(1,2);
% seed outside the image, return
if Xpos < 0.5 ||Ypos < 0.5||Xpos > handles.width||Ypos > handles.height, return, end
% [goodX, goodY] = seedSnapping(hObject, eventdata, handles, Xpos, Ypos, 5);
% Xpos = goodX;
% Ypos = goodY;

tCost = 0; % for debugging window
seedX = Xpos;
seedY = Ypos;
%Left click to put seed
if isempty(handles.XpathPos) %start the first seed point
    handles.XpathPos = Xpos;
    handles.YpathPos = Ypos;
else % A seed point on the way
    [xPath, yPath] = getPath(handles.xPathMap, handles.yPathMap, Xpos, Ypos);
    if isempty(xPath)
        xPath = Xpos;
        yPath = Ypos;
    end
    handles.XpathPos = [handles.XpathPos, double(xPath(:)')];
    handles.YpathPos = [handles.YpathPos, double(yPath(:)')];
end

handles.seedIndx = handles.seedIndx + 1;
handles.seedList(handles.seedIndx) = length(handles.XpathPos); % Save the previous path length for the undo operation
set([handles.contours(handles.maskInd,1), handles.contours(handles.maskInd,2)], ...
    'XData', Xpos, 'YData', Ypos);
drawnow expose
[handles.xPathMap, handles.yPathMap] = determinePath(handles.L, Xpos, Ypos, 200);

% Double-click to close path
if ~(strcmp(get(get(0, 'CurrentFigure'), 'SelectionType'), 'normal'))
    [xPath, yPath] = getPath(handles.xPathMap, handles.yPathMap, handles.XpathPos(1), handles.YpathPos(1));
if isempty(xPath)
    xPath = handles.XpathPos(1);
    yPath = handles.YpathPos(1);
end
handles.XpathPos = [handles.XpathPos, double(xPath(:)')];
handles.YpathPos = [handles.YpathPos, double(yPath(:)')];
set([handles.contours(handles.maskInd,1), handles.contours(handles.maskInd,2)], ...
    'XData', handles.XpathPos, 'YData', handles.YpathPos);
drawnow expose
%mask management
mask = poly2mask(handles.XpathPos, handles.YpathPos,...
    size(handles.current,1), size(handles.current,2));
handles.masks(:,:,handles.maskInd) = mask;
handles.coordinates{handles.maskInd,1} = handles.XpathPos;
handles.coordinates{handles.maskInd,2} = handles.YpathPos;
handles.maskInd = handles.maskInd + 1;

set(get(0, 'CurrentFigure'), 'WindowButtonMotionFcn', '', 'WindowButtonDownFcn', '', 'KeyPressFcn', '');
uiresume(gcf);
end

guidata(hObject, handles);

function mouseMove2(hObject, eventdata, handles)
handles = guidata(hObject);
C = get (gca, 'CurrentPoint');
Xpos = int16(C(1,1));
Ypos = int16(C(1,2));

% no seed point
if (isempty(handles.XpathPos)) || (sum(abs(handles.xPathMap(:))) == 0), return, end
% current position outside of image
if Xpos < 0.5 ||Ypos < 0.5||Xpos > handles.width||Ypos > handles.height, return, end

% Get path from current position to the last seed point
[xPath, yPath] = getPath(handles.xPathMap, handles.yPathMap, Xpos, Ypos);
if isempty(xPath)
    xPath = Xpos;
    yPath = Ypos;
end
set([handles.contours(handles.maskInd,1), handles.contours(handles.maskInd,2)], 'XData', [handles.XpathPos, double(xPath(:)')], ...
                      'YData', [handles.YpathPos, double(yPath(:)')]);
drawnow expose

guidata(hObject, handles);

function keyFuncs1(hObject, eventdata, handles)
handles = guidata(hObject);
switch eventdata.Key
    % ENTER to end path at last seed point
    case 'return'
        [xPath, yPath] = getPath(handles.xPathMap, handles.yPathMap, handles.XpathPos(1), handles.YpathPos(1));
        if isempty(xPath)
            xPath = handles.XpathPos(1);
            yPath = handles.YpathPos(1);
        end
        %handles.XpathPos = [handles.XpathPos, double(xPath(:)')];
        %handles.YpathPos = [handles.YpathPos, double(yPath(:)')];
        set([handles.contours(handles.maskInd,1), handles.contours(handles.maskInd,2)], 'XData', handles.XpathPos, 'YData', handles.YpathPos);
        drawnow expose
        %mask management
        mask = poly2mask(handles.XpathPos, handles.YpathPos,...
            size(handles.current,1), size(handles.current,2));
        handles.masks(:,:,handles.maskInd) = mask;
        handles.coordinates{handles.maskInd,1} = handles.XpathPos;
        handles.coordinates{handles.maskInd,2} = handles.YpathPos;
        handles.maskInd = handles.maskInd + 1;

        set(get(0, 'CurrentFigure'), 'WindowButtonMotionFcn', '', 'WindowButtonDownFcn', '', 'KeyPressFcn', '');
        uiresume(get(0, 'CurrentFigure'));
    % DELETE: delete last seed
    case {'delete', 'backspace'}
        handles.seedIndx = handles.seedIndx - 1;
        if ~handles.seedIndx
            return
        end
        handles.XpathPos = handles.XpathPos(1:handles.seedList(handles.seedIndx));
        handles.YpathPos = handles.YpathPos(1:handles.seedList(handles.seedIndx));
        %delete last length of path, redraw
        set([handles.contours(handles.maskInd,1), handles.contours(handles.maskInd,2)], 'XData', handles.XpathPos, 'YData', handles.YpathPos);
        drawnow expose
        [handles.xPathMap, handles.yPathMap] = determinePath(handles.L, handles.XpathPos(end), ...
            handles.YpathPos(end), 200);
        mouseMove2(hObject, eventdata, handles);
    otherwise
end
guidata(hObject, handles);

% --- Executes on button press in selectContour.
function selectContour_Callback(hObject, eventdata, handles)
% hObject    handle to selectContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curFig = get(0, 'CurrentFigure');
figure(curFig);
try
    set(curFig, 'WindowButtonDownFcn', @getSelectedContour, 'KeyPressFcn', @keyFuncs2,...
        'DoubleBuffer', 'on'); 
catch
    set(curFig, 'WindowButtonDownFcn', @getSelectedContour, 'KeyPressFcn', @keyFuncs2,...
        'DoubleBuffer', 'on');
end
guidata(hObject, handles);

function getSelectedContour(hObject, eventdata, handles)
handles = guidata(hObject);
if size(handles.masks)
    C = get (gca, 'CurrentPoint');
    Xpos = C(1,1);
    Ypos = C(1,2);
    maxInd = handles.maskInd - 1;
    handles.selectedMaskInd = 0;
    for i = 1:maxInd
        if inpolygon(Xpos, Ypos, handles.coordinates{i,1}, handles.coordinates{i,2})
            handles.selectedMaskInd = i;
            if ~isempty(handles.C1)
                delete(handles.C1);
            end
            handles.C1 = plot(handles.coordinates{i,1}, handles.coordinates{i,2}, 'r--','LineWidth',2);
            fprintf('The selected contour index is: %d\n',handles.selectedMaskInd);
            break
        end
    end
    if handles.selectedMaskInd == 0
        fprintf('no contour is selected.\n');
        if ~isempty(handles.C1)
                delete(handles.C1);
        end
    end
else
    h = msgbox('Draw a contour first!');
end

guidata(hObject, handles);

function keyFuncs2(hObject, eventdata, handles)
handles = guidata(hObject);
switch eventdata.Key
    case {'delete', 'backspace'}
        if handles.selectedMaskInd == 0
            h = msgbox('Select a mask first!');
        else
            maxInd = handles.maskInd - 1;
            %delete contour line & rerrange the line array
            delete(handles.contours(handles.selectedMaskInd,:));
            %delete corresponding mask in masks[]
            %delete corresponding coordinates in cell array
            for i = handles.selectedMaskInd:maxInd-1
                handles.masks(:,:,i) = handles.masks(:,:,i+1);
                handles.coordinates(i,:) = handles.coordinates(i+1,:);
                handles.contours(i,:) = handles.contours(i+1,:);
            end
            handles.masks(:,:,maxInd) = [];
            handles.coordinates(maxInd,:) = [];
            %delete highlighted red line
            if ~isempty(handles.C1)
                delete(handles.C1);
            end
            %change maskInd
            handles.selectedMaskInd = 0;
            handles.maskInd = handles.maskInd - 1;
        end
    otherwise
end

guidata(hObject, handles);

% --- Executes on button press in showContours.
function showContours_Callback(hObject, eventdata, handles)
if handles.maskInd - 1 > 0
    for i = 1 : handles.maskInd-1
        set(handles.contours(i,:),'Visible','on');
    end
end
guidata(hObject, handles);

% --- Executes on button press in hideContours.
function hideContours_Callback(hObject, eventdata, handles)
if handles.maskInd - 1 > 0
    for i = 1 : handles.maskInd-1
        set(handles.contours(i,:),'Visible','off');
    end
end
guidata(hObject, handles);


% --- Executes on selection change in debugMode.
function debugMode_Callback(hObject, eventdata, handles)
% hObject    handle to debugMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns debugMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from debugMode

global tCost seedX seedY localCost pq h;

handles = guidata(hObject);
[height, width, dim] = size(handles.current);
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
costGraph = zeros(height*3,width*3,3);
f0 = figure('Name','Debug Window');

switch str{val}
    case 'Pixel Node'
        for i = 1:height
            for j = 1:width
                costGraph(3*(i-1)+2,3*(j-1)+2,:) = handles.current(i,j,:);
            end
        end
        costGraph = uint8(costGraph);      
        image(costGraph);
        axis equal
    case 'Cost Graph'
        localCost = getCost2(handles.current); %local cost, height*width*8 matrix
        for i = 1:height
            for j = 1:width
                % set all center-pixel in each 3*3 window to white
                % costGraph(3*(i-1)+2,3*(j-1)+2,:) = 255; 
                costGraph(3*(i-1)+2,3*(j-1)+2,:) = handles.current(i,j,:);
                costGraph(3*(i-1)+3,3*(j-1)+2,:) = localCost(i,j,1);
                costGraph(3*(i-1)+3,3*(j-1)+1,:) = localCost(i,j,2);
                costGraph(3*(i-1)+2,3*(j-1)+1,:) = localCost(i,j,3);
                costGraph(3*(i-1)+1,3*(j-1)+1,:) = localCost(i,j,4);
                costGraph(3*(i-1)+1,3*(j-1)+2,:) = localCost(i,j,5);
                costGraph(3*(i-1)+1,3*(j-1)+3,:) = localCost(i,j,6);
                costGraph(3*(i-1)+2,3*(j-1)+3,:) = localCost(i,j,7);
                costGraph(3*(i-1)+3,3*(j-1)+3,:) = localCost(i,j,8);
            end
        end
        costGraph = uint8(costGraph);      
        image(costGraph);
        axis equal
        
    case 'Path Tree'
        curX= int16(seedX);
        curY= int16(seedY);      
        localCost = getCost2(handles.current);
        for i = 1:height
            for j = 1:width
                costGraph(3*(i-1)+2,3*(j-1)+2,:) = handles.current(i,j,:);
                costGraph(3*(i-1)+3,3*(j-1)+2,:) = localCost(i,j,1);
                costGraph(3*(i-1)+3,3*(j-1)+1,:) = localCost(i,j,2);
                costGraph(3*(i-1)+2,3*(j-1)+1,:) = localCost(i,j,3);
                costGraph(3*(i-1)+1,3*(j-1)+1,:) = localCost(i,j,4);
                costGraph(3*(i-1)+1,3*(j-1)+2,:) = localCost(i,j,5);
                costGraph(3*(i-1)+1,3*(j-1)+3,:) = localCost(i,j,6);
                costGraph(3*(i-1)+2,3*(j-1)+3,:) = localCost(i,j,7);
                costGraph(3*(i-1)+3,3*(j-1)+3,:) = localCost(i,j,8);
            end
        end
        costGraph = uint8(costGraph);      
        image(costGraph);
        axis equal
        hold on
        plot(3*(curX-1)+2, 3*(curY - 1)+2, 'ro');
        % imageSize
        sizeX = size(localCost,1);
        sizeY = size(localCost,2);
        % expanded/processed state
        state = zeros(sizeX, sizeY); % INITIAL=0, ACTIVE=1, EXPANDED=2
        % path map of x,y
        xPath = zeros(sizeX, sizeY);
        yPath = zeros(sizeX, sizeY);
        % line matrix (yellow and orange ones)
        yline = gobjects(sizeX, sizeY);
        % empty active list (priority queue)
        pq = zeros(3,0); 
        % total cost matrix
        totalCost = Inf(sizeX, sizeY);
        totalCost(curX, curY) = 0.0; %seed with total cost 0
        %insert seed to active list
        insert(curX, curY, 0.0);
        ExpandedCount = 0;
        %%% --- finish initializing, start algorithm --- %%%
        while ~isempty(pq)
            [q_x,q_y,q_cost]=extractmin; % remove minimumcost pixel q from active list
            state(q_x, q_y) = 2; % set as expanded
            ExpandedCount = ExpandedCount +1;                      
            title(gca, ['Expanded: ', num2str(ExpandedCount)]);        
            step = 1;
            for i = 1:8 % loop around neighbors
                [r_x, r_y] = setNeighbor(q_x, q_y, step, i);
                if r_x > 0 && r_y >0 && r_x < sizeX+1 && r_y < sizeY+1
                    if state(r_x,r_y) ~= 2 % r not expanded
                        temptCost = totalCost(q_x,q_y) + localCost(q_x,q_y,i);
                        if (state(r_x, r_y) == 1) && (temptCost < totalCost(r_x, r_y))
                            state(r_x, r_y) = 0;
                        end
                        if state(r_x, r_y) == 0
                            totalCost(r_x, r_y) = temptCost;
                            insert(r_x, r_y,temptCost);
                            state(r_x, r_y) = 1;
                            xPath(r_x, r_y) = q_x;
                            yPath(r_x, r_y) = q_y;                   
                            yline(r_x,r_y) = plot ([(q_x-1)*3+2, (r_x-1)*3+2],[(q_y-1)*3+2, (r_y-1)*3+2], 'y-');
                            drawnow limitrate
                        end
                    end
                end
            end
        end
        % extract overall min, trace back
        while q_x ~= curX || q_y ~= curY
            plot([3*(q_x-1)+2,3*(xPath(q_x,q_y)-1)+2], [3*(q_y-1)+2, 3*(yPath(q_x,q_y)-1)+2],'r-');
            temp_x = xPath(q_x,q_y);
            temp_y = yPath(q_x,q_y);
            q_x = temp_x;
            q_y = temp_y;
            drawnow
        end
        
    case 'Min Path'      
        curX= int16(seedX);
        curY= int16(seedY);
        localCost = getCost2(handles.current);
        for i = 1:height
            for j = 1:width
                costGraph(3*(i-1)+2,3*(j-1)+2,:) = handles.current(i,j,:);
                costGraph(3*(i-1)+3,3*(j-1)+2,:) = localCost(i,j,1);
                costGraph(3*(i-1)+3,3*(j-1)+1,:) = localCost(i,j,2);
                costGraph(3*(i-1)+2,3*(j-1)+1,:) = localCost(i,j,3);
                costGraph(3*(i-1)+1,3*(j-1)+1,:) = localCost(i,j,4);
                costGraph(3*(i-1)+1,3*(j-1)+2,:) = localCost(i,j,5);
                costGraph(3*(i-1)+1,3*(j-1)+3,:) = localCost(i,j,6);
                costGraph(3*(i-1)+2,3*(j-1)+3,:) = localCost(i,j,7);
                costGraph(3*(i-1)+3,3*(j-1)+3,:) = localCost(i,j,8);
            end
        end
        costGraph = uint8(costGraph);      
        image(costGraph);
        axis equal
        hold on
        plot(3*(curX-1)+2, 3*(curY - 1)+2, 'ro');
        h = line('Parent', gca, 'XData', [], 'YData', [], 'Clipping', 'off', ...
    'Color', 'g', 'LineStyle', '-', 'LineWidth', 1.5);
        set(gcf, 'WindowButtonMotionFcn', @mouseMoveDebug);
    otherwise
end

guidata(hObject,handles)

function mouseMoveDebug(hObject, eventdata, handles)
global tCost seedX seedY localCost pq h;
C = get (gca, 'CurrentPoint');
mousePosX = int16(C(1,1));
mousePosY = int16(C(1,2));
curX= int16(seedX);
curY= int16(seedY);

if isPixel(mousePosX, mousePosY)
    %deal with mousePosX&Y to match the original image ones
    mousePosX = (mousePosX - 2)/3 + 1;
    mousePosY = (mousePosY - 2)/3 + 1;
    %updata total cost 
    % imageSize
    sizeX = size(localCost,1);
    sizeY = size(localCost,2);
    % expanded/processed state
    state = zeros(sizeX, sizeY); % INITIAL=0, ACTIVE=1, EXPANDED=2
    % path map of x,y
    xPath = zeros(sizeX, sizeY);
    yPath = zeros(sizeX, sizeY);
    % empty active list (priority queue)
    pq = zeros(3,0); 
    % total cost matrix
    totalCost = Inf(sizeX, sizeY);
    totalCost(curX, curY) = 0.0; %seed with total cost 0
    %insert seed to active list
    insert(curX, curY, 0.0);
    %%% --- finish initializing, start algorithm --- %%%
    while ~isempty(pq) 
        [q_x,q_y,q_cost]=extractmin; % remove minimumcost pixel q from active list
        state(q_x, q_y) = 2; % set as expanded
        %%%%% if meet destination, break
        if q_x == mousePosX && q_y == mousePosY
            break
        end
        step = 1;
        for i = 1:8 % loop around neighbors
            [r_x, r_y] = setNeighbor(q_x, q_y, step, i);
            if r_x > 0 && r_y >0 && r_x < sizeX+1 && r_y < sizeY+1
                if state(r_x,r_y) ~= 2 % r not expanded
                    temptCost = totalCost(q_x,q_y) + localCost(q_x,q_y,i);
                    if (state(r_x, r_y) == 1) && (temptCost < totalCost(r_x, r_y))
                        state(r_x, r_y) = 0;
                    end
                    if state(r_x, r_y) == 0
                        totalCost(r_x, r_y) = temptCost;
                        insert(r_x, r_y,temptCost);
                        state(r_x, r_y) = 1;
                        xPath(r_x, r_y) = q_x;
                        yPath(r_x, r_y) = q_y;
                    end
                end
            end
        end
    end
    %%%%% back trace
    x = [];
    y = [];
    while q_x ~= curX || q_y ~= curY
        x = [x, 3*(q_x-1)+2];
        y = [y, 3*(q_y-1)+2];
        %h = plot([3*(q_x-1)+2,3*(xPath(q_x,q_y)-1)+2], [3*(q_y-1)+2, 3*(yPath(q_x,q_y)-1)+2],'r-');
        set(h,'XData',x,'YData',y);
        temp_x = xPath(q_x,q_y);
        temp_y = yPath(q_x,q_y);
        q_x = temp_x;
        q_y = temp_y;
    end
end

function p = isPixel(i,j)
if mod(i,3) == 2 && mod(j,3) == 2
    p = 1;
else 
    p = 0;
end
% --- Executes during object creation, after setting all properties.
function debugMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to debugMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
