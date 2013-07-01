function varargout = Hill_explorer(varargin)
% HILL_EXPLORER MATLAB code for Hill_explorer.fig
%      HILL_EXPLORER, by itself, creates a new HILL_EXPLORER or raises the existing
%      singleton*.
%
%      H = HILL_EXPLORER returns the handle to a new HILL_EXPLORER or the handle to
%      the existing singleton*.
%
%      HILL_EXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HILL_EXPLORER.M with the given input arguments.
%
%      HILL_EXPLORER('Property','Value',...) creates a new HILL_EXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hill_explorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Hill_explorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Hill_explorer

% Last Modified by GUIDE v2.5 02-Jun-2013 20:24:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Hill_explorer_OpeningFcn, ...
                   'gui_OutputFcn',  @Hill_explorer_OutputFcn, ...
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


% --- Executes just before Hill_explorer is made visible.
function Hill_explorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Hill_explorer (see VARARGIN)

% Choose default command line output for Hill_explorer
handles.output = hObject;

handles.mu=9.537e-4;
handles.ecc=0.048775;
% TODO move these in the draw function
points=2000;
bb=1.5; % Bounding box
x=linspace(-bb,bb,points);
y=linspace(-bb,bb,points);
[xx,yy]=meshgrid(x,y);
handles.xx=xx;
handles.yy=yy;
handles.xlim=[-1.3 1.3];
handles.ylim=[-1.3 1.3];
% handles.xlim=[bb bb];
% handles.ylim=[-bb bb];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Hill_explorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Hill_explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function time_slider_Callback(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
energy=get(handles.energy_slider,'Value');
time=get(hObject,'Value');
set(handles.time_text,'String',num2str(time))
handles.energy=energy;
handles.time=time;
draw(handles);
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function time_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function energy_slider_Callback(hObject, eventdata, handles)
% hObject    handle to energy_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time=get(handles.time_slider,'Value');
energy=get(hObject,'Value');
set(handles.energy_text,'String',num2str(energy))
handles.energy=energy;
handles.time=time;
draw(handles);
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function energy_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energy_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function time_text_Callback(hObject, eventdata, handles)
% hObject    handle to time_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time=str2double(get(hObject,'String'));
set(handles.time_slider,'Value',time);
handles.time=time;
draw(handles)

% Hints: get(hObject,'String') returns contents of time_text as text
%        str2double(get(hObject,'String')) returns contents of time_text as a double


% --- Executes during object creation, after setting all properties.
function time_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function energy_text_Callback(hObject, eventdata, handles)
% hObject    handle to energy_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
energy=str2double(get(hObject,'String'));
set(handles.energy_slider,'Value',energy);
%time=get(handles.time_slider,'Value'); %Fixme
handles.energy=energy;
%handles.time=time; %Fixme
draw(handles)
% Hints: get(hObject,'String') returns contents of energy_text as text
%        str2double(get(hObject,'String')) returns contents of energy_text as a double


% --- Executes during object creation, after setting all properties.
function energy_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energy_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function draw(handles)
h=handles;
zz=Omega(h.xx,h.yy,h.mu)/(1+h.ecc*cos(h.time));
% Code to plot the Potential
%surfc(handles.potential,h.xx,h.yy,-log(zz),'EdgeColor','none')
%title(h.potential,'-log (Omega(x,y,\mu))')
contour(h.hill_region,h.xx,h.yy,zz,[-h.energy,-h.energy],'b');
%title(h.hill_region,'Hill''s region ','fontsize',20)
set(h.hill_region,'fontsize',20)
hz = zoom(h.hill_region);
set(hz,'ActionPostCallback',@mypostcallback);
%set(hz,'Enable','on');
xlim(handles.hill_region,handles.xlim)
ylim(handles.hill_region,handles.ylim)

% --------------------------------------------------------------------
function mypostcallback(obj, evd)
h=guihandles(obj);
old_handles=guidata(h.figure1);
old_handles.xlim=xlim;
old_handles.ylim=ylim;
guidata(h.figure1,old_handles)
