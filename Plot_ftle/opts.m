function varargout = opts(varargin)
%OPTS M-file for opts.fig
%      OPTS, by itself, creates a new OPTS or raises the existing
%      singleton*.
%
%      H = OPTS returns the handle to a new OPTS or the handle to
%      the existing singleton*.
%
%      OPTS('Property','Value',...) creates a new OPTS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to opts_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      OPTS('CALLBACK') and OPTS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in OPTS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help opts

% Last Modified by GUIDE v2.5 02-May-2013 00:06:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @opts_OpeningFcn, ...
                   'gui_OutputFcn',  @opts_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before opts is made visible.
function opts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for opts
flags.lcs=0;
flags.ftle=0;
flags.save_fig=0;
flags.save_mat=0;
flags.spikes=0;
flags.spikes_value=0;
flags.gridfit=0;
flags.gif=0;
handles.flags = flags;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes opts wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = opts_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.flags;
% The figure can be deleted now
delete(handles.figure1);

function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon as text
%        str2double(get(hObject,'String')) returns contents of epsilon as a double


% --- Executes during object creation, after setting all properties.
function epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in spikes_checkbox.
function spikes_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to spikes_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spikes_checkbox
handles.flags.spikes=get(hObject,'Value');
handles.flags.spikes_value=str2double(get(handles.spikes_value,'String'));
guidata(hObject, handles);

function spikes_value_Callback(hObject, eventdata, handles)
% hObject    handle to spikes_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spikes_value as text
%        str2double(get(hObject,'String')) returns contents of spikes_value as a double
set(handles.spikes_checkbox,'Value',1);
handles.flags.spikes=1;
handles.flags.spikes_value=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function spikes_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spikes_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_mat_checkbox.
function save_mat_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_mat_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mat_checkbox
handles.flags.save_mat=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in save_fig_checkbox.
function save_fig_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_fig_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_fig_checkbox
handles.flags.save_fig=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in plot_lcs.
function plot_lcs_Callback(hObject, eventdata, handles)
% hObject    handle to plot_lcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_lcs
handles.flags.lcs=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in gridfit.
function gridfit_Callback(hObject, eventdata, handles)
% hObject    handle to gridfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gridfit
handles.flags.gridfit=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in make_gif.
function make_gif_Callback(hObject, eventdata, handles)
% hObject    handle to make_gif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of make_gif
handles.flags.gif=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close

% --- Executes when user attempts to close figure1.
% FIXME close is different from pushig ok button:
% exiting should mean cancel all, and stop here.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end