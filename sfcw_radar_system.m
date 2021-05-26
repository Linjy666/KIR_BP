function varargout = sfcw_radar_system(varargin)
% SFCW_RADAR_SYSTEM MATLAB code for sfcw_radar_system.fig
%      SFCW_RADAR_SYSTEM, by itself, creates a new SFCW_RADAR_SYSTEM or raises the existing
%      singleton*.
%
%      H = SFCW_RADAR_SYSTEM returns the handle to a new SFCW_RADAR_SYSTEM or the handle to
%      the existing singleton*.
%
%      SFCW_RADAR_SYSTEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SFCW_RADAR_SYSTEM.M with the given input arguments.
%
%      SFCW_RADAR_SYSTEM('Property','Value',...) creates a new SFCW_RADAR_SYSTEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sfcw_radar_system_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sfcw_radar_system_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sfcw_radar_system

% Last Modified by GUIDE v2.5 26-Jun-2017 23:00:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sfcw_radar_system_OpeningFcn, ...
                   'gui_OutputFcn',  @sfcw_radar_system_OutputFcn, ...
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

% --- Executes during object creation, after setting all properties.
function sfcw_radar_system_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sfcw_radar_system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.raw_ad_data = [];
handles.dpd_data = [];
handles.dpd_coef = [];
handles.calibrated_dpd_data = [];
handles.radar_data = [];
handles.syscomp_data = [];
handles.range_profiles = [];
handles.raw_range_profiles = [];
handles.background_range_profiles = [];
handles.declutter_range_profiles = [];
handles.data_info = [];
handles.data_type_text = {'原始中频采样','数字正交解调','参考通道校正',...
    '脉冲压缩','杂波抑制', '脉压校正', '成像', '慢时间维频谱'};
handles.pushbutton_status_array = [];
handles.frame_number = [];
handles.frame_count = [];
handles.ru = [];
handles.flag = [];
handles.movie_flag = 0;
handles.declutter_para = [0,0,40,0.95,0];
handles.declutter_methods = {'two_pulses_canceller',...
    'three_pulses_canceller','sliding_window','exponential_weight','background_subtraction'};
handles.figure_tags = [];
handles.image_info = [];
handles.image_info.axis_mode = 'xoy';
handles.bp_sum = [];
handles.bp_sum_cf = [];
handles.x_tick = [];
handles.y_tick = [];
handles.slow_time_spectrum = [];
handles.slow_time_fft_number = [];
guidata(hObject, handles);

% --- Executes just before sfcw_radar_system is made visible.
function sfcw_radar_system_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sfcw_radar_system (see VARARGIN)

% Choose default command line output for sfcw_radar_system
handles.output = hObject;

% 显示初始化
axes(handles.axes_image);cla;
if ~isempty(handles.data_info)
    close_figures(handles);
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sfcw_radar_system wait for user response (see UIRESUME)
% uiwait(handles.sfcw_radar_system);
% --- Executes during object deletion, before destroying properties.
function sfcw_radar_system_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to sfcw_radar_system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.data_info)
    close_figures(handles);
end


% --- Outputs from this function are returned to the command line.
function varargout = sfcw_radar_system_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% close all;
% clear all;
% clc;
axes(handles.axes_image);cla;
% warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% javaFrame = get(gcf,'javaFrame');
% set(javaFrame,'Maximized',1);

% --------------------------------------------------------------------
function edit_copy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit_file_name,'string'); % 获取文本框中的字符串
clipboard('copy', str); % 将字符串复制到剪切板

% --------------------------------------------------------------------
function edit_past_Callback(hObject, eventdata, handles)
% hObject    handle to edit_past (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste'); % 获取剪切板中的内容
set(handles.edit_file_name,'string',str); % 设置文本框中的内容


% --------------------------------------------------------------------
function edit_cut_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit_file_name,'string'); % 获取文本框中的字符串
clipboard('copy', str); % 将字符串复制到剪切板
set(handles.edit_file_name,'string',''); % 设置文本框中的内容

% --- Executes on button press in pushbutton_open_file.
function pushbutton_open_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
default_path = get(handles.edit_file_name, 'string');
if isempty(default_path)
    default_path = {'*.*'};
else
    [pathstr,name,ext] = fileparts(default_path);
    if isempty(ext)
        default_path = fullfile(default_path, '\*.*');
    end
end
[file_name, path_name] = uigetfile(default_path, '选择数据文件');
if file_name == 0
    return;
end
h = waitbar(0,'处理中...');

[radar_data, data_info] = read_raw_data([path_name, file_name]);

show_data_info([path_name, file_name], data_info);
set(handles.edit_file_name, 'string', [path_name, file_name]);

handles.frame_number = size(radar_data, ndims(radar_data));
handles.data_info = data_info;
handles.ru = 3.0e8/(2*handles.data_info.df);

frame_number = handles.frame_number;
set(handles.slider_frame_count,...
    'max', frame_number, 'min', 1,...
    'SliderStep', [1/(frame_number-1) 10/(frame_number-1)],...
    'value', 1);
set(handles.text_frame_count, 'string', '1');
handles.frame_count = 1;

% 通道成像标签
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number-tr_pair_number;
handles.figure_tags = cell(channel_number,1);
for i=1:ref_channel_number
    handles.figure_tags(i) = cellstr(['Ref Channel', num2str(i)]);
end
if tr_pair_number>1 % 处理多通道的情况
    for i=(1 : tr_pair_number/4)
        handles.figure_tags(i+ref_channel_number) = cellstr(['T1R',num2str(i)]);
    end
    for i=(tr_pair_number/4+1 : 2*tr_pair_number/4)
        handles.figure_tags(i+ref_channel_number) = cellstr(['T2R',num2str(i-tr_pair_number/4)]);
    end
    for i=(2*tr_pair_number/4+1 : 3*tr_pair_number/4)
        handles.figure_tags(i+ref_channel_number) = cellstr(['T3R',num2str(i-2*tr_pair_number/4)]);
    end
    for i=(3*tr_pair_number/4+1 : tr_pair_number)
        handles.figure_tags(i+ref_channel_number) = cellstr(['T4R',num2str(i-3*tr_pair_number/4)]);
    end
end
set(handles.popupmenu_imaging_channel, 'string', handles.figure_tags);

% 数字正交解调系数
sample_number = handles.data_info.sample_number;
dpd_coef = zeros(2, sample_number);
dpd_coef(1, 2:4:end) = 1;
dpd_coef(1, 4:4:end) = -1;
dpd_coef(2, 1:4:end) = -1;
dpd_coef(2, 3:4:end) = 1;
handles.dpd_coef = dpd_coef;

frequency_number = handles.data_info.frequency_number;
switch data_info.data_type
    % 原始中频数据raw
    case 0
        handles.flag = 1;
        fid = fopen('raw_ad.dat', 'wb');
        for i=1:frame_number
            raw_ad_data = reshape(radar_data(:,i), ...
                [sample_number, channel_number, frequency_number]);
            raw_ad_data = permute(raw_ad_data, [1 3 2]);
            fwrite(fid, raw_ad_data, 'double');
        end
        fclose(fid);
        % 调整按钮状态
        handles.pushbutton_status_array = [1,1,2,2,2,2,2,2,2,2];
        adjust_pushbutton_status(handles);
    % 数字正交解调数据dpd
    case 1
        handles.flag = 2;
        fid = fopen('dpd.dat', 'wb');
        fwrite(fid, [real(radar_data), imag(radar_data)], 'double');
        fclose(fid);
        % 调整按钮状态
        handles.pushbutton_status_array = [1,2,1,2,2,2,2,2,2,2];
        adjust_pushbutton_status(handles);
    % 通道校准数据rcc
    case 2
        handles.flag = 3;
        fid = fopen('rcc.dat', 'wb');
        fwrite(fid, [real(radar_data), imag(radar_data)], 'double');
        fclose(fid);
        % 调整按钮状态
        handles.pushbutton_status_array = [1,2,2,1,2,2,2,2,2,2];
        adjust_pushbutton_status(handles);
    % 脉冲压缩数据rp
    case 3
        handles.flag = 4;
        fid = fopen('rp.dat', 'wb');
        fwrite(fid, [real(radar_data), imag(radar_data)], 'double');
        fclose(fid);
        % 调整按钮状态
        handles.pushbutton_status_array = [1,2,2,2,1,2,2,2,2,2];
        adjust_pushbutton_status(handles);
end

show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_dpd.
function pushbutton_dpd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 2;
handles.dpd_data = [];
channel_number = handles.data_info.channel_number;
frequency_number = handles.data_info.frequency_number;
sample_number = handles.data_info.sample_number;

dpd_len = sample_number - mod(sample_number, 4);% 要求为4的整数倍
fid1 = fopen('raw_ad.dat', 'rb');
fid2 = fopen('dpd.dat', 'wb');
h = waitbar(0,'处理中...');
for i=1:handles.frame_number
    raw_ad_data = fread(fid1, [sample_number, frequency_number, channel_number], 'double');
    dpd_data = handles.dpd_coef(:, 1:dpd_len)*reshape(raw_ad_data(1:dpd_len,:,:),dpd_len,[]); % 截取指定dpd_len的中频数据进行dpd
%     dpd_data = complex(dpd_data(2,:), dpd_data(1,:)); % 构建IQ两路复数据
    dpd_data = complex(dpd_data(1,:), dpd_data(2,:)); % 构建IQ两路复数据
    dpd_data = reshape(dpd_data, [frequency_number, channel_number]);
%     handles.dpd_data(:,:,i) = dpd_data;
    fwrite(fid2, [real(dpd_data), imag(dpd_data)], 'double');
    waitbar(i/(handles.frame_number-1),h);
end
fclose(fid1);fclose(fid2);
if strcmpi(get(handles.pushbutton_ref_channel_calibration, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,2,2,2,2,2,2,2];
    adjust_pushbutton_status(handles);
end
show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_ref_channel_calibration.
function pushbutton_ref_channel_calibration_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ref_channel_calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 3;
frame_number = handles.frame_number;
frequency_number = handles.data_info.frequency_number;
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number-tr_pair_number;
H = repmat(hamming(frequency_number), [1, tr_pair_number]);

fid1 = fopen('dpd.dat', 'rb');
fid2 = fopen('rcc.dat', 'wb');
h = waitbar(0,'处理中...');
for i=1:frame_number
    dpd_data = complex(fread(fid1, [frequency_number, channel_number], 'double'),...
        fread(fid1, [frequency_number, channel_number], 'double'));
    ref_channel = dpd_data(:,1:ref_channel_number);
    data_channel = dpd_data(:,ref_channel_number+1:channel_number);
    rcc_data =  data_channel./repmat(ref_channel, [1, tr_pair_number]).*H;
    fwrite(fid2, [real(rcc_data), imag(rcc_data)], 'double');
    waitbar(i/(frame_number-1),h);
end
fclose(fid1);
fclose(fid2);

if strcmpi(get(handles.pushbutton_pulse_compress, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,1,2,2,2,2,2,2];
    adjust_pushbutton_status(handles);
end

show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_pulse_compress.
function pushbutton_pulse_compress_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pulse_compress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 4;
ifft_number = handles.data_info.ifft_number;
frame_number = handles.frame_number;
frequency_number = handles.data_info.frequency_number;
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number-tr_pair_number;

range_profiles = zeros(ifft_number, channel_number);

fid1 = fopen('dpd.dat', 'rb');
fid2 = fopen('rcc.dat', 'rb');
fid3 = fopen('rp.dat', 'wb');
fid4 = fopen('raw_rp.dat', 'wb');
h = waitbar(0,'处理中...');
for i=1:frame_number
    if fid1>=0
        dpd_data = complex(fread(fid1, [frequency_number, channel_number], 'double'),...
            fread(fid1, [frequency_number, channel_number], 'double'));
        range_profiles(:,1:ref_channel_number) = fftshift(ifft(dpd_data(:,1:ref_channel_number), ifft_number));
        raw_range_profiles = ifft(dpd_data, ifft_number);
        fwrite(fid4, [real(raw_range_profiles), imag(raw_range_profiles)], 'double');
    end
    rcc_data = complex(fread(fid2, [frequency_number, channel_number], 'double'),...
        fread(fid2, [frequency_number, channel_number], 'double'));
    range_profiles(:,ref_channel_number+1:channel_number) = ifft(rcc_data, ifft_number);
    fwrite(fid3, [real(range_profiles), imag(range_profiles)], 'double');
    waitbar(i/(frame_number-1),h);
end
if fid1>=0
fclose(fid1);
end
fclose(fid2);fclose(fid3);fclose(fid4);

if strcmpi(get(handles.pushbutton_mitigate_clutter, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,1,1,2,2,2,2,2];
    adjust_pushbutton_status(handles);
end

show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_mitigate_clutter.
function pushbutton_mitigate_clutter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mitigate_clutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 5;
ifft_number = handles.data_info.ifft_number;
frame_number = handles.frame_number;
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number-tr_pair_number;

fid1 = fopen('rp.dat', 'rb');
fid2 = fopen('declutter_rp.dat', 'wb');
h = waitbar(0,'处理中...');
if get(handles.radiobutton_declutter, 'value')
    backgorund = 0;
    index = get(handles.popupmenu_cluter_mitigation, 'value');
    mode = char(handles.declutter_methods(index));
    para = str2double(get(handles.edit_declutter_para, 'string'));
    
    for i=1:frame_number
        range_profiles = complex(fread(fid1, [ifft_number, channel_number], 'double'),...
            fread(fid1, [ifft_number, channel_number], 'double'));
        [declutter_range_profiles, backgorund] = mitigate_clutter(...
            range_profiles(:,ref_channel_number+1:end), backgorund, i, mode, para);
        fwrite(fid2, [real(declutter_range_profiles), imag(declutter_range_profiles)], 'double');
        waitbar(i/(handles.frame_number-1),h);
    end
end

if get(handles.radiobutton_declutter_background, 'value')
    for i=1:frame_number
        range_profiles = complex(fread(fid1, [ifft_number, channel_number], 'double'),...
            fread(fid1, [ifft_number, channel_number], 'double'));
        declutter_range_profiles = range_profiles(:,ref_channel_number+1:end) -...
            handles.background_range_profiles;
        fwrite(fid2, [real(declutter_range_profiles), imag(declutter_range_profiles)], 'double');
        waitbar(i/(frame_number-1),h);
    end
end
fclose(fid1);fclose(fid2);

if strcmpi(get(handles.pushbutton_calibrate_range_profiles, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,1,1,1,2,2,2,2];
    adjust_pushbutton_status(handles);
end

show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_calibrate_range_profiles.
function pushbutton_calibrate_range_profiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calibrate_range_profiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fid = fopen('range_profiles_shift.dat','rb');
if isempty(fid)
    msgbox('脉压校准数据文件不存在！');
    fclose(fid);
    return;
else
    handles.flag = 6;

    frame_number = handles.frame_number;
    ifft_number = handles.data_info.ifft_number;
    tr_pair_number = handles.data_info.tr_pair_number;
    calibrated_range_profiles = zeros(ifft_number, tr_pair_number);
    
%     shift = fread(fid, handles.data_info.tr_pair_number, 'int32');
%     shift = [shift;shift(1:tr_pair_number-length(shift))];
%     shift = [9,14,14,9,7,12,12,7];
    shift = 12*ones(handles.data_info.tr_pair_number,1);
%     fclose(fid);
%     shift = zeros(size(shift));
    fid1 = fopen('declutter_rp.dat', 'rb');
    fid2 = fopen('calibrated_rp.dat', 'wb');
    
    h = waitbar(0,'处理中...');
    handles.calibrated_range_profiles = handles.declutter_range_profiles;
    for i=1:frame_number
        declutter_range_profiles = complex(fread(fid1, [ifft_number, tr_pair_number], 'double'),...
            fread(fid1, [ifft_number, tr_pair_number], 'double'));
        for j=1:tr_pair_number
            calibrated_range_profiles(1:end-shift(j),j) =  declutter_range_profiles(1+shift(j):end,j);
            calibrated_range_profiles(end-shift(j)+1:end,j) = 0;
        end
        fwrite(fid2, [real(calibrated_range_profiles), imag(calibrated_range_profiles)], 'double');
        waitbar(i/(frame_number-1),h);
    end
    fclose(fid1);fclose(fid2);
    show_channel_waveform(handles);
    show_image(handles);
    guidata(hObject,handles);
    close(h);
end
if strcmpi(get(handles.pushbutton_imaging, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,1,1,1,1,1,2,2];
    adjust_pushbutton_status(handles);
end

% --- Executes on button press in pushbutton_imaging.
function pushbutton_imaging_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_imaging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 7;
h = waitbar(0,'处理中...');
if get(handles.radiobutton_xoy_axis, 'value')
    handles.image_info.axis_mode = 'xoy';
end
if get(handles.radiobutton_polar_axis, 'value')
    handles.image_info.axis_mode = 'polar';
end

x_width = str2double(get(handles.popupmenu_x_width, 'string'));
x_width = x_width(get(handles.popupmenu_x_width, 'value'));

dx = str2double(get(handles.popupmenu_dx, 'string'));
dx = dx(get(handles.popupmenu_dx, 'value'));

y_min = str2double(get(handles.popupmenu_y_min, 'string'));
y_min = y_min(get(handles.popupmenu_y_min, 'value'));

y_width = str2double(get(handles.popupmenu_y_width, 'string'));
y_width = y_width(get(handles.popupmenu_y_width, 'value'));

dy = str2double(get(handles.popupmenu_dy, 'string'));
dy = dy(get(handles.popupmenu_dy, 'value'));

z_min = str2double(get(handles.popupmenu_z_min, 'string'));
z_min = z_min(get(handles.popupmenu_z_min, 'value'));

z_width = str2double(get(handles.popupmenu_z_width, 'string'));
z_width = z_width(get(handles.popupmenu_z_width, 'value'));

dz = str2double(get(handles.popupmenu_dz, 'string'));
dz = dz(get(handles.popupmenu_dz, 'value'));

handles.image_info.x_min = -x_width/2;
handles.image_info.dx    = dx;
handles.image_info.x_max = x_width/2;

handles.image_info.y_min = y_min;
handles.image_info.dy    = dy;
handles.image_info.y_max = y_min+y_width;

handles.image_info.z_min = z_min;
handles.image_info.dz    = dz;
handles.image_info.z_max = z_min+z_width;

frame_number = handles.frame_number;
ifft_number = handles.data_info.ifft_number;
tr_pair_number = handles.data_info.tr_pair_number;
dr = handles.ru/ifft_number;
[handles.x_tick, handles.y_tick, handles.z_tick,...
    handles.x_grid, handles.y_grid, handles.z_grid,...
    index_mat, phase_mat] = calculate_3d_image_grid(handles.data_info.fl, dr,...
    handles.data_info.antenna_pos, handles.image_info);

handles.radius = sqrt(handles.x_grid.^2+handles.y_grid.^2+handles.z_grid.^2);
handles.frame_count = round(get(handles.slider_frame_count,'Value'));


fid1 = fopen('calibrated_rp.dat', 'rb');
fid2 = fopen('bp_sum.dat', 'wb');
fid3 = fopen('bp_sum_cf.dat', 'wb');
for i=1:frame_number
    calibrated_range_profiles = complex(fread(fid1, [ifft_number, tr_pair_number], 'double'),...
        fread(fid1, [ifft_number, tr_pair_number], 'double'));
    [bp_sum, bp_sum_cf] = back_projection_3d(...
        calibrated_range_profiles, index_mat, phase_mat);
    fwrite(fid2, bp_sum, 'double');
    fwrite(fid3, bp_sum_cf, 'double');
    waitbar(i/(frame_number-1),h);
end
fclose(fid1);fclose(fid2);fclose(fid3);

if strcmpi(get(handles.pushbutton_movie, 'Enable'), 'off')
    % 调整按钮状态
    handles.pushbutton_status_array = [1,1,1,1,1,1,1,1,1,2];
    adjust_pushbutton_status(handles);
end

show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on button press in pushbutton_slow_time_spectrum.
function pushbutton_slow_time_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_slow_time_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 8;
index = get(handles.popupmenu_slow_time_fft_number, 'value');
mode = get(handles.popupmenu_slow_time_fft_number, 'string');
slow_time_fft_number = str2double(mode(index));
handles.slow_time_fft_number = slow_time_fft_number;
h = waitbar(0,'处理中...');
fid  = fopen('slow_time_spectrum.dat', 'wb');
for i=1:handles.frame_number
    if i<=slow_time_fft_number
        for j=1:handles.data_info.tr_pair_number
            channel_range_profiles = squeeze(handles.calibrated_range_profiles(:,j,1:i));
            handles.slow_time_spectrum = fft(channel_range_profiles, slow_time_fft_number, 2);
            fwrite(fid, handles.slow_time_spectrum, 'double');
        end
    else
        for j=1:handles.data_info.tr_pair_number
            channel_range_profiles = squeeze(handles.calibrated_range_profiles(:,j,i-slow_time_fft_number+1:i));
            handles.slow_time_spectrum = fft(channel_range_profiles, slow_time_fft_number, 2);
            fwrite(fid, handles.slow_time_spectrum, 'double');
        end
    end
    waitbar(i/(handles.frame_number-1),h);
end
fclose(fid);
% for i=1:handles.data_info.tr_pair_number
%     channel_range_profiles = squeeze(handles.calibrated_range_profiles(:,i,:));
%     handles.slow_time_spectrum(:,:,i) = fft(channel_range_profiles, 2048, 2);
%     waitbar(i/(handles.data_info.tr_pair_number-1),h);
% end
handles.radar_data = handles.slow_time_spectrum;
show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);
close(h);

% --- Executes on slider movement.
function slider_frame_count_Callback(hObject, eventdata, handles)
% hObject    handle to slider_frame_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.frame_count = round(get(hObject,'Value'));
set(handles.text_frame_count, 'string', num2str(handles.frame_count));
guidata(hObject,handles);
show_channel_waveform(handles);

data_type_text = char(handles.data_type_text(handles.flag));
switch data_type_text
    case {'成像', '慢时间维频谱'}
%         pushbutton_imaging_Callback(hObject, eventdata, handles);
%     otherwise
        show_image(handles);
end

% --- Executes on selection change in popupmenu_imaging_channel.
function popupmenu_imaging_channel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_imaging_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_imaging_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_imaging_channel
show_image(handles);

% --- Executes on selection change in popupmenu_cluter_mitigation.
function popupmenu_cluter_mitigation_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluter_mitigation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cluter_mitigation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cluter_mitigation
index = get(hObject,'Value');
switch char(handles.declutter_methods(index))
    case 'exponential_weight'
        set(handles.text_declutter_para, 'string', '指数加权系数');
        set(handles.edit_declutter_para, 'string', num2str(handles.declutter_para(index)), 'enable', 'on');
    case 'sliding_window'
        set(handles.text_declutter_para, 'string', '滑窗宽度');
        set(handles.edit_declutter_para, 'string', num2str(handles.declutter_para(index)), 'enable', 'on');
    otherwise
        set(handles.text_declutter_para, 'string', '');
        set(handles.edit_declutter_para, 'string', '', 'enable', 'off');
end

function edit_declutter_para_Callback(hObject, eventdata, handles)
% hObject    handle to edit_declutter_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_declutter_para as text
%        str2double(get(hObject,'String')) returns contents of edit_declutter_para as a double
value = str2double(get(hObject,'String'));
index = get(handles.popupmenu_cluter_mitigation,'Value');
switch char(handles.declutter_methods(index))
    case 'exponential_weight'
        if value>=1 || value <= 0
            msgbox('指数加权系数超出(0, 1)的范围');
            set(hObject, 'string', num2str(handles.declutter_para(index)), 'enable', 'on');
        else
            handles.declutter_para(index) = value;
        end
    case 'sliding_window'
        if value~=fix(value) || value<=0
            msgbox('滑窗宽度必须为正整数');
            set(hObject, 'string', num2str(handles.declutter_para(index)), 'enable', 'on');
        else
            handles.declutter_para(index) = value;
        end
end
guidata(hObject,handles);

function plot_phase_linear_degree(data)
[p, y, m1, n1, m2, n2, min_index, max_index] = evaluate_linear_degree(data);
plot(p);grid on;hold on; plot(y, 'r-');
hold on; plot( [min_index, max_index], [p(min_index), p(max_index)], 'd');
hold on; plot( [min_index, max_index], [y(min_index), y(max_index)], 'rd');
legend('原始相位分布', '拟合相位分布');
annotation(gcf,'textbox',...
    [0.565 0.214 0.44 0.0625],...
    'Color',[1 0 0],...
    'String',...
    sprintf('[%6.2f°, %6.2f°]', wrapTo180(m2), wrapTo180(n2)),...%{['[',num2str(100*m1, '%6.2f'),'% , ', num2str(100*n1, '%6.2f'),'%]']},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');

function show_channel_waveform(handles)
% data = handles.radar_data;
fs = handles.data_info.fs;
sample_number = handles.data_info.sample_number;
frequency_number = handles.data_info.frequency_number;
frame_count = handles.frame_count;
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number - tr_pair_number;
data_type_text = char(handles.data_type_text(handles.flag));
ifft_number = handles.data_info.ifft_number;
figure_tags = handles.figure_tags;
channel_index = get(handles.popupmenu_imaging_channel, 'value');

switch data_type_text
    case {'原始中频采样'}
        fid = fopen('raw_ad.dat','rb');
        fseek(fid, (frame_count-1)*sample_number*frequency_number*channel_number*8,'bof');
        data = fread(fid, [sample_number, frequency_number, channel_number], 'double');
        fclose(fid);
        data = reshape(data, [sample_number*frequency_number, channel_number]);

        show_channel_figures(handles, channel_index);
        subplot(121), plot(data(:,channel_index));grid on;xlabel('样本数据点');
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text]);
        
        spectrum = abs(fft(data(:,channel_index), ifft_number));
        db_spectrum = 20*log10(spectrum/max(spectrum));
        x_tick = (0:ifft_number/2-1)/ifft_number*fs;
        subplot(122), plot(x_tick, db_spectrum(1:ifft_number/2));grid on;xlabel('频率(Hz)');
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text, '频谱']);

    case {'数字正交解调'}
        fid = fopen('dpd.dat','rb');
        fseek(fid, (frame_count-1)*2*frequency_number*channel_number*8,'bof');
        data = complex(fread(fid, [frequency_number, channel_number], 'double'),...
            fread(fid, [frequency_number, channel_number], 'double'));
        fclose(fid);
        show_channel_figures(handles, channel_index);
        subplot(121),plot(abs(data(:,channel_index)));grid on;
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text, '幅度']);
        subplot(122),plot_phase_linear_degree(data(:,channel_index));
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text, '相位']);

    case {'参考通道校正'}
        fid = fopen('rcc.dat','rb');
        fseek(fid, (frame_count-1)*2*frequency_number*tr_pair_number*8,'bof');
        data = complex(fread(fid, [frequency_number, tr_pair_number], 'double'),...
            fread(fid, [frequency_number, tr_pair_number], 'double'));
        fclose(fid);        
        show_channel_figures(handles, channel_index+ref_channel_number);
        subplot(121),plot(abs(data(:,channel_index)));grid on;
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index+ref_channel_number)), data_type_text, '幅度']);
        subplot(122),plot_phase_linear_degree(data(:,channel_index));
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index+ref_channel_number)), data_type_text, '相位']);

    case {'脉冲压缩'}
        fid3 = fopen('dpd.dat','rb');
        fid1 = fopen('rp.dat','rb');
        fid2 = fopen('raw_rp.dat','rb');
        if fid3>=0
            fseek(fid1, (frame_count-1)*2*ifft_number*channel_number*8,'bof');
            data = complex(fread(fid1, [ifft_number, channel_number], 'double'),...
                fread(fid1, [ifft_number, channel_number], 'double'));
            
            fseek(fid2, (frame_count-1)*2*frequency_number*channel_number*8,'bof');
            raw_range_profiles = complex(fread(fid2, [ifft_number, channel_number], 'double'),...
                fread(fid2, [ifft_number, channel_number], 'double'));
            show_channel_figures(handles, channel_index);
            data_abs = abs(data(:,channel_index));
            db_data = 20*log10(data_abs/max(data_abs));
            
            raw_range_profiles_abs = abs(raw_range_profiles(:,channel_index));
            db_raw_range_profiles = 20*log10(raw_range_profiles_abs/max(raw_range_profiles_abs));
            
            x_tick = (0:ifft_number-1)/ifft_number*handles.ru;
            plot(x_tick, db_raw_range_profiles);xlabel('距离(m)');grid on;
            hold on;plot(x_tick, db_data, 'r-');legend('通道校准前', '通道校准后');
            title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text]);
        else
            fseek(fid1, (frame_count-1)*2*ifft_number*channel_number*8,'bof');
            data = complex(fread(fid1, [ifft_number, channel_number], 'double'),...
                fread(fid1, [ifft_number, channel_number], 'double'));
            show_channel_figures(handles, channel_index+ref_channel_number);
            data_abs = abs(data(:,channel_index+ref_channel_number));
            db_data = 20*log10(data_abs/max(data_abs));
            x_tick = (0:ifft_number-1)/ifft_number*handles.ru;
            plot(x_tick, db_data);xlabel('距离(m)');grid on;
            title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index+ref_channel_number)), data_type_text]);
        end
        fclose(fid1);fclose(fid2);
        if (fid3>=0)
            fclose(fid3);
        end
    case {'脉压校正'}
        fid1 = fopen('calibrated_rp.dat','rb');
        fseek(fid1, (frame_count-1)*2*ifft_number*tr_pair_number*8,'bof');
        data = complex(fread(fid1, [ifft_number, tr_pair_number], 'double'),...
            fread(fid1, [ifft_number, tr_pair_number], 'double'));
        show_channel_figures(handles, channel_index);
        data_abs = abs(data(:,channel_index-ref_channel_number));
        db_data = 20*log10(data_abs/max(data_abs));
        x_tick = (0:ifft_number-1)/ifft_number*handles.ru;
        plot(x_tick, db_data);xlabel('距离(m)');grid on;
        title(['第', num2str(frame_count),'帧', char(figure_tags(channel_index)), data_type_text]);
        fclose(fid1);
        
    case {'杂波抑制'}
        fid1 = fopen('declutter_rp.dat','rb');
        fseek(fid1, (frame_count-1)*2*ifft_number*tr_pair_number*8,'bof');
        data = complex(fread(fid1, [ifft_number, tr_pair_number], 'double'),...
            fread(fid1, [ifft_number, tr_pair_number], 'double'));
        show_channel_figures(handles, channel_index);
        data_abs = abs(data(:,channel_index-ref_channel_number));
        db_data = 20*log10(data_abs/max(data_abs));
        x_tick = (0:ifft_number-1)/ifft_number*handles.ru;
        plot(x_tick, db_data);xlabel('距离(m)');grid on;
        title(['第', num2str(frame_count),'帧', char(handles.figure_tags(channel_index)), data_type_text]);
        fclose(fid1);
%     case {'慢时间维频谱'}
%         fid  = fopen('slow_time_spectrum.dat', 'rb');
%         frame_count = round(get(handles.slider_frame_count,'Value'));
%         frame_size = handles.data_info.ifft_number * handles.slow_time_fft_number * 8;
%         tr_pair_number = handles.data_info.tr_pair_number;
%         slow_time_fft_number = handles.slow_time_fft_number;
%         position = (frame_count-1)*tr_pair_number*frame_size;
%         for i=ref_channel_number+1:handles.data_info.channel_number
%             show_channel_figures(handles, i);
%             fseek(fid, (i-ref_channel_number-1)*frame_size + position, 'bof');
%             data = fread(fid, [handles.data_info.ifft_number,handles.slow_time_fft_number], 'double'); 
%             
%             data_abs = abs(data);
%             db_data = 20*log10(data_abs(:,1:slow_time_fft_number/2)/max(max(data_abs(:,1:slow_time_fft_number/2))));
%             
%             slow_time_fft_number = size(data_abs, 2);
%             x_tick = (0:slow_time_fft_number/2-1)*handles.data_info.prf/slow_time_fft_number;
%             y_tick=(0:handles.data_info.ifft_number-1)/handles.data_info.ifft_number*handles.ru;
%             imagesc(x_tick, y_tick, db_data);grid on;colormap('jet');
%             xlabel('频率(Hz)');ylabel('距离(m)');colorbar;
%             title(['第', num2str(frame_count),'帧', char(handles.figure_tags(i)), data_type_text]);
%         end
%         fclose(fid);
end

function show_image(handles)
axes(handles.axes_image);
frame_number = handles.frame_number;
frequency_number = handles.data_info.frequency_number;
sample_number = handles.data_info.sample_number;
ifft_number = handles.data_info.ifft_number;
channel_number = handles.data_info.channel_number;
tr_pair_number = handles.data_info.tr_pair_number;
ref_channel_number = channel_number - tr_pair_number;
data_type_text = char(handles.data_type_text(handles.flag));
figure_tags = handles.figure_tags;
channel_index = get(handles.popupmenu_imaging_channel, 'value');

%'原始中频采样','数字正交解调','参考通道校正','脉冲压缩','杂波抑制','脉压校正','成像'
switch data_type_text
    case {'原始中频采样'}
        data = zeros(sample_number*frequency_number, frame_number);
        fid = fopen('raw_ad.dat','rb');
        for i=1:frame_number
            raw_ad = fread(fid, [sample_number, frequency_number, channel_number], 'double');            
            raw_ad = reshape(raw_ad, [sample_number*frequency_number, channel_number]);
            data(:,i) = raw_ad(:,channel_index);
        end
        fclose(fid);
        
        imagesc(abs(data)); title([char(figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);

    case {'数字正交解调'}
        data = zeros(frequency_number, frame_number);
        fid = fopen('dpd.dat','rb');
        for i=1:frame_number
            dpd = complex(fread(fid, [frequency_number, channel_number], 'double'),...
                fread(fid, [frequency_number, channel_number], 'double'));
            data(:,i) = dpd(:,channel_index);
        end
        fclose(fid);

        imagesc(abs(data)); title([char(figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);
        
    case {'参考通道校正'}
        if channel_index <= ref_channel_number
            return;
        else
            data = zeros(frequency_number, frame_number);
            fid = fopen('rcc.dat','rb');
            for i=1:frame_number
                rcc = complex(fread(fid, [frequency_number, channel_number], 'double'),...
                    fread(fid, [frequency_number, channel_number], 'double'));
                data(:,i) = rcc(:,channel_index-ref_channel_number);
            end
            fclose(fid);
            imagesc(abs(data)); title([char(figure_tags(channel_index-ref_channel_number)), data_type_text]);
        end
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);
        
    case {'脉冲压缩'}
        if channel_index <= ref_channel_number
            if~isempty(handles.dpd_data)
                data = zeros(ifft_number, frame_number);
                fid = fopen('rp.dat','rb');
                for i=1:frame_number
                    rp = complex(fread(fid, [ifft_number, channel_number], 'double'),...
                        fread(fid, [ifft_number, channel_number], 'double'));
                    data(:,i) = rp(:,channel_index);
                end
                fclose(fid);
%                 data = squeeze(handles.range_profiles(:,channel_index,:));
            else
                return;
            end
        else
            data = zeros(ifft_number, frame_number);
            fid = fopen('rp.dat','rb');
            for i=1:frame_number
                rp = complex(fread(fid, [ifft_number, tr_pair_number], 'double'),...
                    fread(fid, [ifft_number, tr_pair_number], 'double'));
                data(:,i) = rp(:,channel_index);
            end
            fclose(fid);
        end
        x_tick = 1:size(data,2);y_tick=(0:ifft_number-1)/ifft_number*handles.ru;
        imagesc(x_tick, y_tick, abs(data));ylabel('距离(m)');
        title([char(figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);

    case {'脉压校正'}
        data = zeros(ifft_number, frame_number);
        fid = fopen('calibrated_rp.dat','rb');
        for i=1:frame_number
            rp = complex(fread(fid, [ifft_number, tr_pair_number], 'double'),...
                fread(fid, [ifft_number, tr_pair_number], 'double'));
            data(:,i) = rp(:,channel_index);
        end
        fclose(fid);
        if channel_index <= ref_channel_number
            data = squeeze(handles.range_profiles(:,channel_index,:));
        else
            data = data(:,channel_index-ref_channel_number);
        end
        x_tick = 1:size(data,2);y_tick=(0:ifft_number-1)/ifft_number*handles.ru;
        imagesc(x_tick, y_tick, abs(data));ylabel('距离(m)');
        title([char(handles.figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);
        
    case {'杂波抑制'}
        data = zeros(ifft_number, frame_number);
        fid = fopen('declutter_rp.dat','rb');
        for i=1:frame_number
            rp = complex(fread(fid, [ifft_number, channel_number], 'double'),...
                fread(fid, [ifft_number, channel_number], 'double'));
            data(:,i) = rp(:,channel_index);
        end
        fclose(fid);
        if channel_index <= ref_channel_number
            data = squeeze(handles.range_profiles(:,channel_index,:));
        else
            data = data(:,channel_index-ref_channel_number);
        end
        x_tick = 1:size(data,2);y_tick=(0:ifft_number-1)/ifft_number*handles.ru;
        imagesc(x_tick, y_tick, abs(data));ylabel('距离(m)');
        title([char(figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('探测次数');colorbar;
        show_channel_waveform(handles);
        
    case {'成像'}
        frame_count = round(get(handles.slider_frame_count,'Value'));
        fid = fopen('bp_sum_cf.dat','rb');
        [x, y, z] = size(handles.x_grid);
        fseek(fid, (frame_count-1)*x*y*z*8, 'bof');
        data = reshape(fread(fid, x*y*z, 'double'), [x,y,z]);
        fclose(fid);
        show_3d_image(handles.x_grid, handles.y_grid, handles.z_grid, abs(data));
        xlim([min(handles.x_tick) max(handles.x_tick)]);
        ylim([min(handles.y_tick) max(handles.y_tick)]);
        zlim([min(handles.z_tick) max(handles.z_tick)]);
        grid on;
        if strcmpi(handles.image_info.axis_mode,'xoy')
            xlabel('方位(米)');ylabel('距离(米)');zlabel('俯仰(米)')
        end
        if strcmpi(handles.image_info.axis_mode,'polar')
            xlabel('方位(度)');ylabel('距离(米)');zlabel('俯仰(米)')
        end
        title(['第',num2str(handles.frame_count),'帧 ',...
            handles.image_info.axis_mode,'成像结果']);
    case {'慢时间维频谱'}
        fid  = fopen('slow_time_spectrum.dat', 'rb');
        channel_index = get(handles.popupmenu_imaging_channel, 'value');
        if channel_index <= ref_channel_number
            return;
%         else
%             data = squeeze(handles.radar_data(:,:,channel_index-ref_channel_number));
        end
        frame_count = round(get(handles.slider_frame_count,'Value'));
        frame_size = handles.data_info.ifft_number * handles.slow_time_fft_number * 8;
        tr_pair_number = handles.data_info.tr_pair_number;
        slow_time_fft_number = handles.slow_time_fft_number;
        position = (frame_count-1)*tr_pair_number*frame_size;
        fseek(fid, (channel_index-ref_channel_number-1)*frame_size + position, 'bof');
        data = fread(fid, [handles.data_info.ifft_number,handles.slow_time_fft_number], 'double');
        fclose(fid);
%         slow_time_fft_number = size(handles.radar_data, 2);
%         x_tick = (0:slow_time_fft_number/2-1)*handles.data_info.prf/slow_time_fft_number;
%         y_tick=(0:handles.data_info.ifft_number-1)/handles.data_info.ifft_number*handles.ru;
%         imagesc(x_tick, y_tick, abs(data(:,1:slow_time_fft_number/2)));
        x_tick = (0:slow_time_fft_number-1)*handles.data_info.prf/slow_time_fft_number;
        y_tick=(0:handles.data_info.ifft_number-1)/handles.data_info.ifft_number*handles.ru;
        imagesc(x_tick, y_tick, abs(data));
        title([char(handles.figure_tags(channel_index)), data_type_text]);
        grid on;colormap('jet'); xlabel('频率(Hz)');ylabel('距离(m)');colorbar;        
end

function show_channel_figures(handles, i)
h = findobj('tag', char(handles.figure_tags(i)));
if ~isempty(h)
    figure(h); clf;
else
    figure('tag', char(handles.figure_tags(i)));
end

function close_figures(handles)
for i=1:handles.data_info.channel_number
    h = findobj('tag', char(handles.figure_tags(i)));
    if ~isempty(h)
        close(h);
    end
end
h = findobj('tag', 'data_info_figure');
if ~isempty(h)
    close(h);
end

function adjust_pushbutton_status(handles)
array = handles.pushbutton_status_array;
status = {'on', 'off', 'inactive'};
set(handles.pushbutton_show_file_data,          'Enable', char(status(array(1))));
set(handles.pushbutton_dpd,                     'Enable', char(status(array(2))));
set(handles.pushbutton_ref_channel_calibration, 'Enable', char(status(array(3))));
set(handles.pushbutton_pulse_compress,          'Enable', char(status(array(4))));
set(handles.pushbutton_mitigate_clutter,        'Enable', char(status(array(5))));
set(handles.pushbutton_calibrate_range_profiles,'Enable', char(status(array(6))));
set(handles.pushbutton_imaging, 'Enable', char(status(array(7))));
set(handles.pushbutton_slow_time_spectrum, 'Enable', char(status(array(8))));
set(handles.pushbutton_movie,   'Enable', char(status(array(9))));
set(handles.pushbutton_pause,   'Enable', char(status(array(10))));
switch handles.data_info.data_type
    % 数字正交解调数据dpd
    case 1
        % 调整按钮状态
        set(handles.pushbutton_dpd,                     'Enable', 'off');
    % 通道校准数据rcc
    case 2
        % 调整按钮状态
        set(handles.pushbutton_dpd,                     'Enable', 'off');
        set(handles.pushbutton_ref_channel_calibration, 'Enable', 'off');
    % 脉冲压缩数据rp
    case 3
        % 调整按钮状态
        set(handles.pushbutton_dpd,                     'Enable', 'off');
        set(handles.pushbutton_ref_channel_calibration, 'Enable', 'off');
        set(handles.pushbutton_pulse_compress,          'Enable', 'off');
end

% --- Executes on button press in radiobutton_declutter.
function radiobutton_declutter_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_declutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_declutter
set(handles.pushbutton_background_data, 'enable', 'off');
set(handles.popupmenu_cluter_mitigation, 'enable', 'on');
index = get(handles.popupmenu_cluter_mitigation, 'Value');
switch char(handles.declutter_methods(index))
    case {'exponential_weight','sliding_window'}
        set(handles.edit_declutter_para, 'enable', 'on');
end

% --- Executes on button press in radiobutton_declutter_background.
function radiobutton_declutter_background_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_declutter_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_declutter_background
set(handles.popupmenu_cluter_mitigation, 'enable', 'off');
set(handles.pushbutton_background_data, 'enable', 'on');
set(handles.edit_declutter_para, 'enable', 'off');


% --- Executes on button press in pushbutton_background_data.
function pushbutton_background_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_background_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
default_path = get(handles.edit_file_name, 'string');
if isempty(default_path)
    default_path = {'*.*'};
else
    [pathstr,name,ext] = fileparts(default_path);
    if isempty(ext)
        default_path = fullfile(default_path, '\*.*');
    end
end
[file_name, path_name] = uigetfile(default_path, '选择数据文件');
% [file_name, path_name] = uigetfile({'*.*'},'选择背景数据文件');
if file_name == 0
    return;
end
h = waitbar(0,'处理中...');
[radar_data, data_info] = read_raw_data([path_name, file_name]);
frame_number = size(radar_data, 3);
radar_data = radar_data(:,:,fix(frame_number/2));
channel_number = data_info.channel_number;
tr_pair_number = data_info.tr_pair_number;
ref_channel_number = channel_number-tr_pair_number;
frequency_number = data_info.frequency_number;
sample_number = data_info.sample_number;
ifft_number = data_info.ifft_number;

handles.background_range_profiles = [];
% 提取数据
switch data_info.data_type
    % 原始中频数据
    case 0
        radar_data = reshape(radar_data, ...
            [sample_number, channel_number, frequency_number]);
        raw_ad_data = permute(radar_data, [1 3 2]);
        % 数字正交解调
        
        dpd_len = data_info.sample_number;
        dpd_coef = zeros(2, dpd_len);
        dpd_coef(1, 2:4:end) = 1;
        dpd_coef(1, 4:4:end) = -1;
        dpd_coef(2, 1:4:end) = -1;
        dpd_coef(2, 3:4:end) = 1;
        H = repmat(hamming(frequency_number), [1, tr_pair_number]);

        dpd_data = dpd_coef*reshape(raw_ad_data(:,:,:), sample_number,[]); % 截取指定dpd_len的中频数据进行dpd
        % dpd_data = complex(dpd_data(2,:), dpd_data(1,:)); % 构建IQ两路复数据
        dpd_data = complex(dpd_data(1,:), dpd_data(2,:)); % 构建IQ两路复数据
        dpd_data = reshape(dpd_data, [frequency_number, channel_number]);
        
        % 通道校正
        ref_channel = dpd_data(:,1:ref_channel_number);
        data_channel = dpd_data(:,ref_channel_number+1:channel_number);
        rcc_data =  data_channel./repmat(ref_channel, [1, tr_pair_number]).*H;
        
        % 脉冲压缩
        handles.background_range_profiles = ifft(rcc_data, ifft_number);
            
    % 数字正交解调数据
    case 1
        dpd_data = radar_data;
        H = repmat(hamming(frequency_number), [1, tr_pair_number]);
        % 通道校正
        ref_channel = dpd_data(:,1:ref_channel_number);
        data_channel = dpd_data(:,ref_channel_number+1:channel_number);
        rcc_data =  data_channel./repmat(ref_channel, [1, tr_pair_number]).*H;
        
        % 脉冲压缩
        handles.background_range_profiles = ifft(rcc_data, ifft_number);
            
    % 通道校准数据RCC
    case 2
        rcc_data = radar_data;
        handles.background_range_profiles = ifft(rcc_data, ifft_number);% 脉冲压缩        
    % 脉冲压缩数据
    case 3
        handles.background_range_profiles = radar_data;
end

guidata(hObject,handles);
close(h);


% --- Executes on button press in radiobutton_xoy_axis.
function radiobutton_xoy_axis_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_xoy_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_xoy_axis
handles.image_info.axis_mode = 'xoy';
set(handles.popupmenu_x_width, 'string', num2cell(10:2:30));
set(handles.popupmenu_x_width, 'value', 4);
set(handles.popupmenu_dx, 'string', num2cell([0.01,0.02,0.04,0.05,0.08,0.1]));
set(handles.popupmenu_dx, 'value', 2);
set(handles.popupmenu_angle, 'visible', 'on');
set(handles.text_angle, 'visible', 'on');
guidata(hObject,handles);

% --- Executes on button press in radiobutton_polar_axis.
function radiobutton_polar_axis_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_polar_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_polar_axis
handles.image_info.axis_mode = 'polar';
set(handles.popupmenu_x_width, 'string', num2cell(80:20:180));
set(handles.popupmenu_x_width, 'value', 3);
set(handles.popupmenu_dx, 'string', num2cell([0.1,0.2,0.4,0.5,0.8,1.0]));
set(handles.popupmenu_dx, 'value', 2);
set(handles.popupmenu_angle, 'visible', 'off');
set(handles.text_angle, 'visible', 'off');
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_file_data.
function pushbutton_show_file_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show_file_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.data_info.data_type
    % 原始中频数据raw
    case 0
        handles.flag = 1;
        handles.radar_data = reshape(handles.raw_ad_data,...
            [handles.data_info.sample_number*handles.data_info.frequency_number,...
            handles.data_info.channel_number, handles.frame_number]);
    % 数字正交解调数据dpd
    case 1
        handles.flag = 2;
        handles.radar_data = handles.dpd_data;
    % 通道校准数据rcc
    case 2
        handles.flag = 3;
        handles.radar_data = handles.calibrated_dpd_data;
    % 脉冲压缩数据rp
    case 3
        handles.flag = 4;
        handles.radar_data = handles.range_profiles;
end

show_channel_waveform(handles);
show_image(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_movie.
function pushbutton_movie_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 7;
% 调整按钮状态
handles.pushbutton_status_array = [1,1,1,1,1,1,1,1,1,1];
adjust_pushbutton_status(handles);
global movie_flag;
movie_flag = 1;
handles.movie_flag = 1;
fps = str2double(get(handles.edit_fps, 'string'));
frame_count = round(get(handles.slider_frame_count,'Value'));
for i=frame_count:handles.frame_number
    set(handles.slider_frame_count,'Value', i);
    set(handles.text_frame_count,'string', num2str(i));
    show_image(handles);
    pause(1/fps);
    if movie_flag == 0
        break;
    end
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton_pause.
function pushbutton_pause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global movie_flag;
movie_flag = 0;
guidata(hObject,handles);



% --- Executes on button press in pushbutton_save_figure.
function pushbutton_save_figure_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_f_handle = figure('visible', 'off');
new_axes = copyobj(handles.axes_image, new_f_handle); %picture是GUI界面绘图的坐标系句柄
set(new_axes,'units','default','position','default');
[filename, pathname, fileindex] = uiputfile({'*.fig';'*.png';'*.jpg';'*.bmp'}, '保存图像为');
if ~filename
     return
else
      file = fullfile(pathname, filename);
      switch fileindex %根据不同的选择保存为不同的类型
      case 1
          set(new_f_handle, 'visible', 'on'); 
      end
end
saveas(new_f_handle, file);
delete(new_f_handle);
