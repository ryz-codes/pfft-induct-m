function varargout = test_tables_gui(varargin)
% TEST_TABLES_GUI MATLAB code for test_tables_gui.fig
%      TEST_TABLES_GUI, by itself, creates a new TEST_TABLES_GUI or raises the existing
%      singleton*.
%
%      H = TEST_TABLES_GUI returns the handle to a new TEST_TABLES_GUI or the handle to
%      the existing singleton*.
%
%      TEST_TABLES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_TABLES_GUI.M with the given input arguments.
%
%      TEST_TABLES_GUI('Property','Value',...) creates a new TEST_TABLES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_tables_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_tables_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_tables_gui

% Last Modified by GUIDE v2.5 09-Oct-2012 13:30:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_tables_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @test_tables_gui_OutputFcn, ...
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

function update(handles)
% extract the end points
z_min = handles.g.z_min
z_rng = handles.g.z_max - z_min

% make range
d=z_min+z_rng*get(handles.sldD,'Value');
z=z_min+z_rng*get(handles.sldZ,'Value');
r=0:0.001:0.1;

s = sprintf('z=%g, zp=%g',z,d);
disp(s)

[Gr] = handles.g.lookup(r,z,d);
[Gfdm r2] = analy_sol_hs(z+d,handles.L);
[G3 r3] = fdm_run(z,d,handles.L);

cla(handles.axes1);
hold(handles.axes1,'on');
plot(handles.axes1,r,real(Gr),'or')
plot(handles.axes1,r3,real(G3),'--r')
plot(handles.axes1,r2,real(Gfdm))
hold(handles.axes1,'off');
legend(handles.axes1,'Reconstructed solution','FDM solution','Analytical solution');
title(handles.axes1,['Re(G) ' s]);
xlim(handles.axes1,[0,0.1])

cla(handles.axes2);
hold(handles.axes2,'on');
plot(handles.axes2,r,imag(Gr),'or')
plot(handles.axes2,r3,imag(G3),'--r')
plot(handles.axes2,r2,imag(Gfdm))
hold(handles.axes2,'off');
legend(handles.axes2,'Reconstructed solution','FDM solution','Analytical solution');
title(handles.axes2,['Im(G) ' s]);
xlim(handles.axes2,[0,0.1])


% --- Executes just before test_tables_gui is made visible.
function test_tables_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_tables_gui (see VARARGIN)

% Choose default command line output for test_tables_gui
handles.output = hObject;

% Load our pregenerated TPH object
load('default1.mat');
handles.g = g;
handles.L = g.L;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_tables_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_tables_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sldZ_Callback(hObject, eventdata, handles)
% hObject    handle to sldZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sldZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sldZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sldD_Callback(hObject, eventdata, handles)
% hObject    handle to sldD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sldD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sldD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pshUPDATE.
function pshUPDATE_Callback(hObject, eventdata, handles)
% hObject    handle to pshUPDATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update(handles);
