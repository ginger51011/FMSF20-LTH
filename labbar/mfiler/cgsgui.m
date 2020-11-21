function varargout = cgsgui(varargin)
%cgsgui Grafiskt anvandargranssnitt som illustrerar centrala gransvardessatsen
%
% cgsgui
%
% Grafiskt anvandargranssnitt for centrala gransvardessatsen. Illusterar
% hur summa och medelvarde av ett antal standardfordelningar beter sig.
%
% Exempel:
%   cgsgui

% Copyright 2017 Johan Lindström, Olof Zetterqvist, Lena Zetterqvist; Lunds Universitet.

% Last Modified by GUIDE v2.5 05-Sep-2015 18:46:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @cgsgui_OpeningFcn, ...
  'gui_OutputFcn',  @cgsgui_OutputFcn, ...
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


% --- Executes just before cgsgui is made visible.
function cgsgui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cgsgui (see VARARGIN)

handles.fordelning = 'Fördelning';
handles.mu = 0.4;
handles.sigma = 0.6;
handles.N = 1;
handles.normplot = 0;
handles.medel = 0;

set(handles.checkbox3,'value',1);

set(handles.edit1,'Enable','off');
set(handles.edit2,'Enable','off');
set(handles.edit1,'String',handles.mu);
set(handles.edit2,'String',handles.sigma);
set(handles.edit3,'String',handles.N);

axes(handles.axes1);
title('Sannolikhets/t\"athetsfunktion', 'interpreter', 'latex');
axes(handles.axes2);
title('Histogram \"over 1000 simulerade summor/medelv\"arden', ...
  'interpreter', 'latex');

% Choose default command line output for cgsgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cgsgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cgsgui_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fordelning.
function fordelning_Callback(hObject, ~, handles)
% hObject    handle to fordelning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fordelning contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fordelning
contents = cellstr(get(hObject,'String'));
current = contents(get(hObject,'Value'));
handles.fordelning = current;

if strcmp(current, 'Fördelning')
  set(handles.text1,'string','Parameter 1');
  set(handles.text2,'string','Parameter 2');
  set(handles.edit1,'Enable','off');
  set(handles.edit2,'Enable','off');
elseif strcmp(current, 'Exp(1/mu)')
  set(handles.text1,'string','mu');
  set(handles.text2,'string','Parameter 2');
  set(handles.edit1,'Enable','on');
  set(handles.edit2,'Enable','off');
  
  handles.mu = 1;
  handles.sigma = 0;
  set(handles.edit1,'String',handles.mu);
  set(handles.edit2,'String',handles.sigma);
  set(handles.edit3,'String',handles.N);
elseif strcmp(current, 'Re(a,b)')
  set(handles.text1,'string','a');
  set(handles.text2,'string','b');
  set(handles.edit1,'Enable','on');
  set(handles.edit2,'Enable','on');
  
  handles.mu = 0;
  handles.sigma = 1;
  set(handles.edit1,'String',handles.mu);
  set(handles.edit2,'String',handles.sigma);
  set(handles.edit3,'String',handles.N);
elseif strcmp(current,'N(mu,sigma)')
  set(handles.text1,'string','mu');
  set(handles.text2,'string','sigma');
  set(handles.edit1,'Enable','on');
  set(handles.edit2,'Enable','on');
  
  handles.mu = 0;
  handles.sigma = 1;
  set(handles.edit1,'String',handles.mu);
  set(handles.edit2,'String',handles.sigma);
  set(handles.edit3,'String',handles.N);
elseif strcmp(current, 'Bin(n,p)')
  set(handles.text1,'string','n');
  set(handles.text2,'string','p');
  set(handles.edit1,'Enable','on');
  set(handles.edit2,'Enable','on');
  
  handles.mu = 10;
  handles.sigma = 0.5;
  set(handles.edit1,'String',handles.mu);
  set(handles.edit2,'String',handles.sigma);
  set(handles.edit3,'String',handles.N);
elseif strcmp(current, 'Po(mu)')
  set(handles.text1,'string','mu');
  set(handles.text2,'string','Parameter 2');
  set(handles.edit1,'Enable','on');
  set(handles.edit2,'Enable','off');
  
  handles.mu = 1;
  handles.sigma = 0;
  set(handles.edit1,'String',handles.mu);
  set(handles.edit2,'String',handles.sigma);
  set(handles.edit3,'String',handles.N);
end
if ~strcmp(handles.fordelning,'Fördelning')
  handles = set_data(handles);
  handles = plot_figures(handles);
  guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function fordelning_CreateFcn(hObject, ~, ~)
% hObject    handle to fordelning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(~, ~, ~)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over fordelning.
function fordelning_ButtonDownFcn(~, ~, ~)
% hObject    handle to fordelning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, ~)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%---------------------------------------------------------------------

function handles = plot_figures(handles)
axes(handles.axes1);
if strcmp(handles.fordelning, 'Bin(n,p)') || strcmp(handles.fordelning, 'Po(mu)')
  bar(handles.x,handles.y,0.2);
elseif ~strcmp(handles.fordelning,'Fördelning')
  plot(handles.x,handles.y);
end
hold on;
plot(handles.varden,zeros(1,length(handles.varden)),'xr');
plot(mean(handles.varden),0,'o');
hold off;


if ((strcmp(handles.fordelning,'Bin(n,p)') && handles.N*handles.mu*handles.sigma*(1-handles.sigma) < 10) || (strcmp(handles.fordelning,'Po(mu)') && handles.mu*handles.N < 5))
  if strcmp(handles.fordelning,'Bin(n,p)')
    temp = binoinv(0.001,handles.mu*handles.N,handles.sigma):binoinv(0.999,handles.mu*handles.N,handles.sigma);
  else
    temp = poissinv(0.001,handles.mu*handles.N):poissinv(0.999,handles.mu*handles.N);
  end
else
  temp = 15;
end

axes(handles.axes2);
if handles.medel
  if handles.N ~= 1
    hist(mean(handles.varden2),temp);
    [handles.varde, handles.x2] = hist(mean(handles.varden2),temp);
  else
    hist(handles.varden2,temp);
    [handles.varde, handles.x2] = hist(handles.varden2,temp);
  end
else
  if handles.N ~= 1
    hist(sum(handles.varden2),temp);
    [handles.varde, handles.x2] = hist(sum(handles.varden2),temp);
  else
    hist(handles.varden2,temp);
    [handles.varde, handles.x2] = hist(handles.varden2,temp);
  end
end

current = handles.fordelning;
if strcmp(current, 'Exp(1/mu)')
  if handles.medel
    mu = handles.mu;
    sigma = handles.mu / sqrt(handles.N);
  else
    mu = handles.mu*handles.N;
    sigma = handles.mu * sqrt(handles.N);
  end
elseif strcmp(current, 'Re(a,b)')
  if handles.medel
    mu = (handles.mu + handles.sigma)/2;
    sigma = (handles.sigma - handles.mu) /sqrt(handles.N)/sqrt(12);
  else
    mu = (handles.mu + handles.sigma)*handles.N/2;
    sigma = (handles.sigma - handles.mu) *sqrt(handles.N/12);
  end
elseif strcmp(current,'N(mu,sigma)')
  if handles.medel
    mu = handles.mu;
    sigma = handles.sigma/sqrt(handles.N);
  else
    mu = handles.mu*handles.N;
    sigma = handles.sigma*sqrt(handles.N);
  end
elseif strcmp(current, 'Bin(n,p)')
  if handles.medel
    mu = handles.mu*handles.sigma;
    sigma = sqrt(handles.mu*handles.sigma*(1 - handles.sigma)/handles.N);
  else
    mu = handles.mu*handles.sigma*handles.N;
    sigma = sqrt(handles.N)*sqrt(handles.mu*handles.sigma*(1 - handles.sigma));
  end
elseif strcmp(current, 'Po(mu)')
  if handles.medel
    mu = handles.mu;
    sigma = sqrt(handles.mu/handles.N);
  else
    mu = handles.mu*handles.N;
    sigma = sqrt(handles.mu*handles.N);
  end
end

axes(handles.axes1);
title('Sannolikhets/t\"athetsfunktion', 'interpreter', 'latex');
axes(handles.axes2);
title('Histogram \"over 1000 simulerade summor/medelv\"arden', ...
  'interpreter', 'latex');
xlabel(sprintf('E = %2.2f, D = %2.2f', mu, sigma), 'interpreter', 'latex');

%Normalanpassning
if handles.normplot
  set(handles.axes3,'Visible','on');
  %figure 1
  x = linspace(mu - 3*sigma,mu + 3*sigma,1000);
  y = normpdf(x,mu,sigma)*sum((handles.x2(2)-handles.x2(1)).*handles.varde);
  axes(handles.axes2);
  hold on;
  plot(x,y,'r','LineWidth',2);
  hold off;
  
  %figure 2
  axes(handles.axes3);
  if handles.medel
    normplot(mean(handles.varden2,1));
  else
    normplot(sum(handles.varden2,1));
  end
  title('Normalf\"ordelning', 'interpreter', 'latex')
  xlabel(''); ylabel('');
  p = [.001 .01 0.05 0.25 .5 0.75 0.95 .99 .999];
  set(handles.axes3, 'yticklabel', p, 'ytick', norminv(p));
else
  cla(handles.axes3, 'reset')
  set(handles.axes3, 'Visible', 'off');
end

function handles = set_data(handles)
current = handles.fordelning;
if strcmp(current, 'Exp(1/mu)')
  handles.varden = exprnd(handles.mu,handles.N,1);
  handles.varden2 = exprnd(handles.mu,handles.N,1000);
  handles.x = linspace(0,expinv(0.999,handles.mu),1000);
  handles.y = exppdf(handles.x,handles.mu);
elseif strcmp(current, 'Re(a,b)')
  handles.varden = unifrnd(handles.mu,handles.sigma,handles.N,1);
  handles.varden2 = unifrnd(handles.mu,handles.sigma,handles.N,1000);
  handles.x = linspace(handles.mu-2,handles.sigma+2,1000);
  handles.y = unifpdf(handles.x,handles.mu,handles.sigma);
elseif strcmp(current,'N(mu,sigma)')
  
  handles.varden = normrnd(handles.mu,handles.sigma,handles.N,1);
  handles.varden2 = normrnd(handles.mu,handles.sigma,handles.N,1000);
  handles.x = linspace(handles.mu - 4*handles.sigma,handles.mu + 4*handles.sigma,1000);
  handles.y = normpdf(handles.x,handles.mu,handles.sigma);
elseif strcmp(current, 'Bin(n,p)')
  handles.varden = binornd(handles.mu,handles.sigma,handles.N,1);
  handles.varden2 = binornd(handles.mu,handles.sigma,handles.N,1000);
  handles.x = 0:handles.mu+2;
  handles.y = binopdf(handles.x,handles.mu,handles.sigma);
elseif strcmp(current, 'Po(mu)')
  handles.x = 0:poissinv(0.999,handles.mu);
  handles.y = poisspdf(handles.x,handles.mu);
  handles.varden = poissrnd(handles.mu,handles.N,1);
  handles.varden2 = poissrnd(handles.mu,handles.N,1000);
end



function edit1_Callback(hObject, ~, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
if ~strcmp(handles.fordelning,'Fördelning')
  handles.mu = str2double(get(hObject,'String'));
  handles = set_data(handles);
  handles = plot_figures(handles);
  guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, ~, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
if ~strcmp(handles.fordelning,'Fördelning')
  handles.sigma = str2double(get(hObject,'String'));
  handles = set_data(handles);
  handles = plot_figures(handles);
  guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, ~, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
if ~strcmp(handles.fordelning,'Fördelning')
  handles.N = floor(str2double(get(hObject,'String')));
  if handles.N < 1
    handles.N = 1;
    set(handles.edit3,'String',1);
  end
  handles = set_data(handles);
  handles = plot_figures(handles);
  guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton.
function pushbutton_Callback(~, ~, handles)
% hObject    handle to pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current = handles.fordelning;
if strcmp(current, 'Exp(1/mu)')
  handles.varden2 = exprnd(handles.mu,1,1000);
elseif strcmp(current, 'Re(a,b)')
  handles.varden2 =unifrnd(handles.mu,handles.sigma,1,1000);
elseif strcmp(current,'N(mu,sigma)')
  handles.varden2 =normrnd(handles.mu,handles.sigma,1,1000);
elseif strcmp(current, 'Bin(n,p)')
  handles.varden2 = binornd(handles.mu,handles.sigma,1,1000);
elseif strcmp(current, 'Po(mu)')
  handles.varden2 =poissrnd(handles.mu,1,1000);
end
for n = 2:handles.N
  handles.N = n;
  if strcmp(current, 'Exp(1/mu)')
    handles.varden2 = [handles.varden2; exprnd(handles.mu,1,1000)];
  elseif strcmp(current, 'Re(a,b)')
    handles.varden2 = [handles.varden2; unifrnd(handles.mu,handles.sigma,1,1000)];
  elseif strcmp(current,'N(mu,sigma)')
    handles.varden2 = [handles.varden2 ; normrnd(handles.mu,handles.sigma,1,1000)];
  elseif strcmp(current, 'Bin(n,p)')
    handles.varden2 = [handles.varden2 ; binornd(handles.mu,handles.sigma,1,1000)];
  elseif strcmp(current, 'Po(mu)')
    handles.varden2 = [handles.varden2 ; poissrnd(handles.mu,1,1000)];
  end
  handles = plot_figures(handles);
  set(handles.edit3,'String',n);
  pause(0.5);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~strcmp(handles.fordelning,'Fördelning')
  handles = set_data(handles);
  handles = plot_figures(handles);
  guidata(hObject, handles);
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, ~, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

handles.normplot = get(hObject,'Value');
if ~strcmp(handles.fordelning,'Fördelning')
  handles = plot_figures(handles);
end
guidata(hObject, handles);



% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, ~, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

handles.medel = 1;
set(handles.checkbox3,'value',0);
if ~strcmp(handles.fordelning,'Fördelning')
  handles = plot_figures(handles);
end
guidata(hObject, handles);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, ~, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
handles.medel = 0;
set(handles.checkbox2,'value',0);
if ~strcmp(handles.fordelning,'Fördelning')
  handles = plot_figures(handles);
end
guidata(hObject, handles);
