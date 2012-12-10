function varargout = run_me_gui(varargin)
% RUN_ME_GUI M-file for run_me_gui.fig
%      RUN_ME_GUI, by itself, creates a new RUN_ME_GUI or raises the existing
%      singleton*.
%
%      H = RUN_ME_GUI returns the handle to a new RUN_ME_GUI or the handle to
%      the existing singleton*.
%
%      RUN_ME_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUN_ME_GUI.M with the given input arguments.
%
%      RUN_ME_GUI('Property','Value',...) creates a new RUN_ME_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before run_me_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to run_me_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% TO DOS
% - Shifting lower bounds. Changing lower bounds changes the offsets.
% - Implement add layer, delete layer, move up and move down.
% - Design the T+H splitting into the gui.

% Edit the above text to modify the response to help run_me_gui

% Last Modified by GUIDE v2.5 08-Oct-2012 21:46:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @run_me_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @run_me_gui_OutputFcn, ...
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

%--------------------------------------------------------------------------
% Refresh functions
%--------------------------------------------------------------------------
% --- Refreshes the display to match the state variables
function disp_refresh(handles)
% handles    handle to figure
L = handles.L;
ii = handles.lay_select;

% Update frequency
set(handles.txtW,'String',num2str(L.w/2/pi));

% Update the list box
listtxt = cell(L.layerN,1);
for ij = 1:L.layerN
    lay_text = L.names{L.layerN-ij+1};
    if isempty(lay_text)
        lay_text = sprintf('Layer %u ',ij);
    end
    listtxt{ij} = lay_text;
end

set(handles.listLAYER,'Value',L.layerN-ii+1)
set(handles.listLAYER,'String',listtxt)

% Update the layer info panel
set(handles.txtNAME,'String',L.names{ii})
set(handles.txtMU_R,'String',num2str(L.mu_r(ii),6))
set(handles.txtSIG,'String',num2str(L.sig(ii),6))
set(handles.txtZN,'String',num2str(L.zN(ii),6))
set(handles.txtBNDS,'String',num2str(L.bnds(ii),6))
set(handles.chkCOIL,'Value',(L.coil_layer == ii))
set(handles.listUNITS,'Value',handles.type(ii)+1)

% Thickness
set(handles.txtBND2,'String',handles.thick(ii))

% Draw rectangles onto plot representing layers
cla(handles.axes1);
for ij = 1:L.layerN
    x = -1;
    y = L.bnds(ij);
    w = 2;
    h = L.bnds(ij+1)-L.bnds(ij);
    
    if L.mu_r(ij) == 1 && L.sig(ij) == 0
    % draw as a white rectangle with black border if it is air
        rectangle('parent',handles.axes1,'position',[x y w h], ...
            'FaceColor',[0.95 0.95 0.95]);
    else
    % draw as grey rectangles without borders if it is not air
        if mod(ij,2)==0
            rectangle('parent',handles.axes1,'position',[x y w h], ...
                'FaceColor',[0.8 0.8 0.8],'LineStyle','none');
        else
            rectangle('parent',handles.axes1,'position',[x y w h], ...
                'FaceColor',[0.7 0.7 0.7],'LineStyle','none');
        end
    end
end
% Draw gridlines if it is turned on
if get(handles.tgGRID,'Value')
    for ij = 1:L.layerN
        z = linspace(L.bnds(ij),L.bnds(ij+1),L.zN(ij));
        y = kron(z(2:end-1).',[1 1]);
        x = kron(ones(size(y,1),1),[-1 1]);
        line(x.',y.','parent',handles.axes1,'color',[0.5,0.5,0.5],'LineWidth',0.1);
        
        bny = [z(1) z(end); z(1) z(end)];
        bnx = x(1:2,1:2).';
        
        line(bnx,bny,'parent',handles.axes1,'color',[0,0,0],'LineWidth',0.6);
    end
end
% Draw text to annotate the rectangles, if labels are turned on
if get(handles.tgLABELS,'Value')
    for ij = 1:L.layerN
        y = L.bnds(ij);
        h = L.bnds(ij+1)-L.bnds(ij);
        cx = 0*(-1)^mod(ij,2);
        cy = y+h/2;

        s = L.names{ij};
        
        % Add details if they are turned on
        if get(handles.tgDETAILS,'Value') 
            if L.coil_layer == ij
                s = [s sprintf('\n (Coil layer)')];
            end

            s = [s sprintf('\n\\mu_r = %3g, \\sigma = %3g, N_z = %3g',...
                L.mu_r(ij),L.sig(ij),L.zN(ij))];
        else
            if L.coil_layer == ij
                s = [s '*'];
            end
        end
        
        text(cx,cy,s,'HorizontalAlignment','center');
    end
end
xlim(handles.axes1,[-1 1]);

pad = 0.1*(L.bnds(end)-L.bnds(1)); % 10% padding
ylim(handles.axes1,[L.bnds(1)-pad L.bnds(end)+pad])


% --- Refreshes the state variables to match data entered
function handles = data_refresh(hObject,handles)

% handles handle to figure
L = handles.L;
ii = handles.lay_select;

% Update frequency
L.w = 2*pi*str2double(get(handles.txtW,'String'));

% Update the layer info panel
L.names{ii} = get(handles.txtNAME,'String');
L.mu_r(ii) = str2double(get(handles.txtMU_R,'String'));
L.sig(ii) = str2double(get(handles.txtSIG,'String'));
L.zN(ii) = str2double(get(handles.txtZN,'String'));

% Change coil layer only if it can be changed
if get(handles.chkCOIL,'Value')
    L.coil_layer = ii;
end
handles.L = L;

% Boundary processing
this_bnd = str2double(get(handles.txtBNDS,'String'));
handles.thick(ii) = str2double(get(handles.txtBND2,'String'));
handles = update_bnds(handles,ii,this_bnd);

% Skin depth processing
handles = refresh_skin_d(handles);

% Save data
guidata(hObject,handles);


function handles = refresh_skin_d(handles)
ii = handles.lay_select;
handles.skin_d = sqrt(2./(handles.L.w .* handles.L.mu_r .* handles.L.mu .* handles.L.sig));
new_type = get(handles.listUNITS,'Value')-1;
if new_type ~= handles.type(ii) % Do we require unit conversion?
    if isfinite(handles.skin_d(ii))
        if new_type == 0
            handles.thick(ii) = handles.thick(ii)*handles.skin_d(ii);
        else
            handles.thick(ii) = handles.thick(ii)/handles.skin_d(ii);
        end
        handles.type(ii) = new_type;    
    else
        waitfor(warndlg('Warning: This layer does not have a finite skin depth.'));
        set(handles.listUNITS,'Value',1);
    end
end



%--------------------------------------------------------------------------
% Layer manipulation
%--------------------------------------------------------------------------
function handles = layer_new(handles,ij)
% Makes new layer by inserting at the end and copying everything in ii.
L = handles.L;

% New layer info
ii = L.layerN+1;
L.layerN = ii;
L.mu_r(ii) = L.mu_r(ij);
L.sig(ii) = L.sig(ij);
L.zN(ii) = L.zN(ij);
L.names{ii} = ['Copy of ' L.names{ij}];
handles.L = L;

% Skin depths
handles.type(ii) = handles.type(ij);
handles = refresh_skin_d(handles);

% Treat bounds
handles.thick(ii) = handles.thick(ij);
handles=update_bnds(handles,ij,L.bnds(ij));

function handles = layer_del(handles,ii)

if handles.L.layerN == 1
    waitfor(warndlg('You only have one layer left! Cannot delete the last layer'));
    return
end

Ln = handles.L;
L =  handles.L;
Ln.mu_r=L.mu_r([1:(ii-1) (ii+1):end]);
Ln.sig=L.sig([1:(ii-1) (ii+1):end]);
Ln.zN=L.zN([1:(ii-1) (ii+1):end]);
Ln.names={L.names{1:(ii-1)} L.names{(ii+1):end}};
Ln.layerN = L.layerN-1;
handles.L = Ln;




if ii==1 % Delete layer 1

    % Skin depths
    handles.type = handles.type([1:(ii-1) (ii+1):end]);
    handles = refresh_skin_d(handles);

    % treat bounds
    handles.thick = handles.thick([1:(ii-1) (ii+1):end]);

    handles=update_bnds(handles,ii,L.bnds(ii));
    
else
    % Adjust the layer you are looking at.
    handles.lay_select = handles.lay_select -1;

    % Skin depths
    handles.type = handles.type([1:(ii-1) (ii+1):end]);
    handles = refresh_skin_d(handles);

    % treat bounds
    handles.thick = handles.thick([1:(ii-1) (ii+1):end]);
    handles=update_bnds(handles,ii-1,L.bnds(ii-1));
    
    if L.coil_layer == ii
        L.coil_layer = ii-1;
    end
end

% Treat coil layer
if ii == L.coil_layer
    Ln.coil_layer=L.coil_layer-1;
end


function handles = shift_swap(handles,orig,dest)
L = handles.L;
Ln = handles.L;

% New layer info
Ln.mu_r(dest) = L.mu_r(orig);
Ln.sig(dest) = L.sig(orig);
Ln.zN(dest) = L.zN(orig);
Ln.names{dest} = L.names{orig};
Ln.mu_r(orig) = L.mu_r(dest);
Ln.sig(orig) = L.sig(dest);
Ln.zN(orig) = L.zN(dest);
Ln.names{orig} = L.names{dest};

handles.L = Ln;

% Skin depths & Layer thickness
temp1 = handles.type(dest); temp2 = handles.thick(dest);
handles.type(dest) = handles.type(orig);
handles.thick(dest) = handles.thick(orig);
handles.type(orig) = temp1;
handles.thick(orig) = temp2;

handles = refresh_skin_d(handles);
handles=update_bnds(handles,1,L.bnds(1));

function handles = update_bnds(handles,ii,this_bnd)
% Updates the bnds properties to match those state variables stored in bnds
thick = handles.thick;

new_bnd = zeros(1,length(thick)+1);
for ij = 1:length(thick)
    if handles.type(ij) == 1
        new_bnd(ij+1) = new_bnd(ij) + thick(ij)*handles.skin_d(ij);
    else
        new_bnd(ij+1) = new_bnd(ij) + thick(ij);
    end
end
dbnd = this_bnd-new_bnd(ii);
new_bnd = new_bnd+dbnd;

handles.L.bnds = new_bnd;


%--------------------------------------------------------------------------
% Toolbars, Buttons and the Select List
%--------------------------------------------------------------------------
function uisave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uisave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = uiputfile('layers.mat','Save layer structure variable as');
if FileName ~= 0 % Didn't cancel
    L = handles.L;
    save(FileName,'L');
end

function uiopen_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uiopen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = uigetfile('layers.mat','Open layer structure variable');
if FileName ~= 0 % Didn't cancel
    load(FileName);
    handles.L = L;
    % TO DO: FOOL PROOF THE LOADING PROCESS
end
handles.lay_select = 1;
handles.thick = handles.L.bnds(2:end)-handles.L.bnds(1:end-1);
handles.type = zeros(handles.L.layerN,1);
handles = refresh_skin_d(handles);

handles.lay_select = 1;
disp_refresh(handles);
guidata(hObject, handles);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = uiputfile('tables.mat','Save Greens function table as');
if FileName ~= 0 % Didn't cancel
    ubnd = str2double(inputdlg('Chose a table upper bound','dialog',1,{'0.99'}));
    lbnd = str2double(inputdlg('Chose a table lower bound','dialog',1,{'0.01'}));
    if (~isnan(ubnd) && ~isnan(lbnd) && ubnd < 1 && lbnd < 1 && ubnd > 0 ...
            && lbnd > 0 && ubnd > lbnd)
        g = tph(handles.L,lbnd,ubnd);
        save(FileName,'g');
    end 
end


% --- Executes on button press in pshADD.
function pshADD_Callback(hObject, eventdata, handles)
% hObject    handle to pshADD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = handles.lay_select;
handles = layer_new(handles,ii);
handles.lay_select=handles.L.layerN;
% refresh displays
disp_refresh(handles);
guidata(hObject,handles);


% --- Executes on button press in pshMINUS.
function pshMINUS_Callback(hObject, eventdata, handles)
% hObject    handle to pshMINUS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = handles.lay_select;
handles = layer_del(handles,ii);
if ii > handles.L.layerN, handles.lay_select=handles.L.layerN; end
% refresh displays
disp_refresh(handles);
guidata(hObject,handles);

% --- Executes on button press in pshUP.
function pshUP_Callback(hObject, eventdata, handles)
% hObject    handle to pshUP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = handles.lay_select;
if ii<handles.L.layerN
    handles=shift_swap(handles,ii,ii+1);
    handles.lay_select=ii+1;
    disp_refresh(handles);
    guidata(hObject,handles);
end


% --- Executes on button press in pshDOWN.
function pshDOWN_Callback(hObject, eventdata, handles)
% hObject    handle to pshDOWN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = handles.lay_select;
if ii > 1
    handles = shift_swap(handles,ii,ii-1);
    handles.lay_select=ii-1;
    disp_refresh(handles);
    guidata(hObject,handles);
end

% --- Executes on selection change in listLAYER.
function listLAYER_Callback(hObject, eventdata, handles)
% hObject    handle to listLAYER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = data_refresh(hObject,handles);
rev_lay = get(hObject,'Value');
handles.lay_select = handles.L.layerN-rev_lay+1;
disp_refresh(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Form Interface
%--------------------------------------------------------------------------
% --- Executes just before run_me_gui is made visible.
function run_me_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_me_gui (see VARARGIN)

% Choose default command line output for run_me_gui
handles.output = hObject;

% Initialize data
handles.L = defaultL(4);
handles.lay_select = 1;
handles.thick = handles.L.bnds(2:end)-handles.L.bnds(1:end-1);
handles.type = zeros(handles.L.layerN,1);
L = handles.L;
handles = refresh_skin_d(handles);

% Refresh 
disp_refresh(handles);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = run_me_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--------------------------------------------------------------------------
% Data entry Callback and Creation
%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function listLAYER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listLAYER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMU_R_Callback(hObject, eventdata, handles)
% hObject    handle to txtMU_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtMU_R as text
%        str2double(get(hObject,'String')) returns contents of txtMU_R as a double


% --- Executes during object creation, after setting all properties.
function txtMU_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMU_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSIG_Callback(hObject, eventdata, handles)
% hObject    handle to txtSIG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtSIG as text
%        str2double(get(hObject,'String')) returns contents of txtSIG as a double


% --- Executes during object creation, after setting all properties.
function txtSIG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSIG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZN_Callback(hObject, eventdata, handles)
% hObject    handle to txtZN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtZN as text
%        str2double(get(hObject,'String')) returns contents of txtZN as a double


% --- Executes during object creation, after setting all properties.
function txtZN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtBNDS_Callback(hObject, eventdata, handles)
% hObject    handle to txtBNDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtBNDS as text
%        str2double(get(hObject,'String')) returns contents of txtBNDS as a double


% --- Executes during object creation, after setting all properties.
function txtBNDS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBNDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtBND2_Callback(hObject, eventdata, handles)
% hObject    handle to txtBND2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtBND2 as text
%        str2double(get(hObject,'String')) returns contents of txtBND2 as a double


% --- Executes during object creation, after setting all properties.
function txtBND2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBND2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkCOIL.
function chkCOIL_Callback(hObject, eventdata, handles)
% hObject    handle to chkCOIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% refresh displays
handles = data_refresh(hObject,handles);
disp_refresh(handles);

% Hint: get(hObject,'Value') returns toggle state of chkCOIL





function txtNAME_Callback(hObject, eventdata, handles)
% hObject    handle to txtNAME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);

% Hints: get(hObject,'String') returns contents of txtNAME as text
%        str2double(get(hObject,'String')) returns contents of txtNAME as a double

% --- Executes during object creation, after setting all properties.
function txtNAME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNAME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in tgGRID.
function tgGRID_Callback(hObject, eventdata, handles)
% hObject    handle to tgGRID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hint: get(hObject,'Value') returns toggle state of tgGRID


% --- Executes on button press in tgDETAILS.
function tgDETAILS_Callback(hObject, eventdata, handles)
% hObject    handle to tgDETAILS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hint: get(hObject,'Value') returns toggle state of tgDETAILS


% --- Executes on button press in tgLABELS.
function tgLABELS_Callback(hObject, eventdata, handles)
% hObject    handle to tgLABELS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hint: get(hObject,'Value') returns toggle state of tgLABELS


% --- Executes on selection change in listUNITS.
function listUNITS_Callback(hObject, eventdata, handles)
% hObject    handle to listUNITS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);


% Hints: contents = cellstr(get(hObject,'String')) returns listUNITS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listUNITS


% --- Executes during object creation, after setting all properties.
function listUNITS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listUNITS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtW_Callback(hObject, eventdata, handles)
% hObject    handle to txtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = data_refresh(hObject,handles);
% refresh displays
disp_refresh(handles);
% Hints: get(hObject,'String') returns contents of txtW as text
%        str2double(get(hObject,'String')) returns contents of txtW as a double


% --- Executes during object creation, after setting all properties.
function txtW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
