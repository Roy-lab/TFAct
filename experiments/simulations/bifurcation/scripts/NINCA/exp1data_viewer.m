function varargout = exp1data_viewer(varargin)
% EXP1DATA_VIEWER M-file for exp1data_viewer.fig
%      EXP1DATA_VIEWER, by itself, creates a new EXP1DATA_VIEWER or raises the existing
%      singleton*.
%
%      H = EXP1DATA_VIEWER returns the handle to a new EXP1DATA_VIEWER or the handle to
%      the existing singleton*.
%
%      EXP1DATA_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXP1DATA_VIEWER.M with the given input arguments.
%
%      EXP1DATA_VIEWER('Property','Value',...) creates a new EXP1DATA_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exp1data_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exp1data_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exp1data_viewer

% Last Modified by GUIDE v2.5 29-Mar-2011 00:29:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exp1data_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @exp1data_viewer_OutputFcn, ...
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

% --- Executes just before exp1data_viewer is made visible.
function exp1data_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exp1data_viewer (see VARARGIN)

% Choose default command line output for exp1data_viewer
handles.output = hObject;

% Process passed in arguments
for i = 1:2:length(varargin)-1
    switch ( varargin{i} )
        case 'K'
            handles.K = varargin{i+1};
        case 'S'
            handles.S = varargin{i+1};
        case 'tfa_index'
            handles.tfa_index = varargin{i+1};
        otherwise
            warning('argument ''%s'' unrecognized',varargin{i});
    end
end

% Update handles structure
guidata(hObject, handles);


% This sets up the initial plot - only do when we are invisible
% so window can get raised using exp1data_viewer.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes exp1data_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exp1data_viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

tp = get(handles.popupmenu1, 'Value');
tfa = get(handles.popupmenu2, 'Value');

% read tfa into buffer ( convert 1x1xK -> 1xK )
S_buffer = zeros(4,handles.K); % exp x K
for exp = 1:4
    for k = 1:handles.K
        if ( handles.tfa_index(exp,tfa) ~= 0 )
            S_buffer(exp,k) = handles.S(exp,handles.tfa_index(exp,tfa),k);
        end
    end
end

hold on;

switch (tp)
    case 1 % elutriation
        t = 0:30:390;
        if ( handles.tfa_index(1,tfa) ~= 0 )
            plot( t, S_buffer(1,1:14), 'b-');
        end
        if ( handles.tfa_index(2,tfa) ~= 0 )
            plot( t, S_buffer(2,1:14), 'b:');
        end
        if ( handles.tfa_index(3,tfa) ~= 0 )
            plot( t, S_buffer(3,1:14), 'b-.');
        end
        if ( handles.tfa_index(4,tfa) ~= 0 )
            plot( t, S_buffer(4,1:14), 'b--');
        end
        
        % set axis nicely
%         v = axis;
%         axis([0 400 v(3) v(4)]);
        axis tight;
        
    case 2 % \alpha-factor
        t = 0:7:119;
        if ( handles.tfa_index(1,tfa) ~= 0 )
            plot( t, S_buffer(1,15:32), 'b-');
        end
        if ( handles.tfa_index(2,tfa) ~= 0 )
            plot( t, S_buffer(2,15:32), 'b:');
        end
        if ( handles.tfa_index(3,tfa) ~= 0 )
            plot( t, S_buffer(3,15:32), 'b-.');
        end
        if ( handles.tfa_index(4,tfa) ~= 0 )
            plot( t, S_buffer(4,15:32), 'b--');
        end
        
        % set axis nicely
%         v = axis;
%         axis([0 120 v(3) v(4)]);
        axis tight;
        
    case 3 % cdc-15
        t = [10 30 50 70:10:250 270 290];
        if ( handles.tfa_index(1,tfa) ~= 0 )
            plot( t, S_buffer(1,33:56), 'b-');
        end
        if ( handles.tfa_index(2,tfa) ~= 0 )
            plot( t, S_buffer(2,33:56), 'b:');
        end
        if ( handles.tfa_index(3,tfa) ~= 0 )
            plot( t, S_buffer(3,33:56), 'b-.');
        end
        if ( handles.tfa_index(4,tfa) ~= 0 )
            plot( t, S_buffer(4,33:56), 'b--');
        end
        
        % set axis nicely
%         v = axis;
%         axis([0 300 v(3) v(4)]);
        axis tight;
        
    case 4 % \alpha-factor (mutant)
        t = 0:15:180;
        if ( handles.tfa_index(1,tfa) ~= 0 )
            plot( t, S_buffer(1,57:69), 'b-');
        end
        if ( handles.tfa_index(2,tfa) ~= 0 )
            plot( t, S_buffer(2,57:69), 'b-');
        end
        if ( handles.tfa_index(3,tfa) ~= 0 )
            plot( t, S_buffer(3,57:69), 'b-');
        end
        if ( handles.tfa_index(4,tfa) ~= 0 )
            plot( t, S_buffer(4,57:69), 'b-');
        end
        
        % set axis nicely
%         v = axis;
%         axis([0 180 v(3) v(4)]);
        axis tight;
        
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
set(hObject, 'String', {'1 elutriation','2 alpha-factor','3 cdc-15'});



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',{'1 Ace2','2 Fkh1','3 Fkh2','4 Mbp1','5 Mcm1',...
                      '6 Ndd1','7 Skn7','8 Stb1','9 Swi4','10 Swi5','11 Swi6'});


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
