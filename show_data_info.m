function show_data_info(file_name, data_info)

waveform_name = {'����Ƶ������', '�������'};
data_type_name = {'ԭʼ��Ƶ����','�����������','ͨ��У׼����','����ѹ������'};
str_data_info(1)  = cellstr(['������ַ��     ', data_info.ip_addr(1:15)]);
str_data_info(2)  = cellstr(['�źŲ��Σ�     ', char(waveform_name(data_info.waveform+1))]);
str_data_info(3)  = cellstr(['�������ͣ�     ', char(data_type_name(data_info.data_type+1))]);
str_data_info(4)  = cellstr(['��ʼƵ�ʣ�     ', num2str(data_info.fl/1e9), ' GHz']);
str_data_info(5)  = cellstr(['��ֹƵ�ʣ�     ', num2str(data_info.fh/1e9), ' GHz']);
str_data_info(6)  = cellstr(['Ƶ�ʲ�����     ', num2str(data_info.df/1e6), ' MHz']);
str_data_info(7)  = cellstr(['������Ƶ��     ', num2str(data_info.prf), ' Hz']);
str_data_info(8)  = cellstr(['��Ƶ������     ', num2str(data_info.fs/1e6), ' MHz']);
str_data_info(9)  = cellstr(['Ƶ��������     ', num2str(data_info.frequency_number)]);
str_data_info(10)  = cellstr(['Ƶ������������ ', num2str(data_info.sample_number)]);
str_data_info(11) = cellstr(['ȫ��ͨ�������� ', num2str(data_info.channel_number)]);
str_data_info(12) = cellstr(['�շ�ͨ�������� ', num2str(data_info.tr_pair_number)]);
str_data_info(13) = cellstr(['���ջ���  �棺 ', num2str(data_info.receiver_gain)]);
str_data_info(14) = cellstr(['IFFT�任������ ', num2str(data_info.ifft_number)]);

h = findobj('tag', 'data_info_figure');
if ~isempty(h),
    figure(h); clf;
else
    h = figure('name', '�״������ļ���Ϣ',...
        'tag','data_info_figure', ...
        'numbertitle','off', ...
        'position',[500   80   500   500], ...
        'menubar', 'none');
end
data_info_text = '�״������ļ�ͷ';
title_string = file_name;
c1 = 0.5;    c2 = 0.5;    c3 = 0.5;
whitebg(h, [c1 c2 c3]);
set(gca,'visible','off');%,'SortMethod','fast');
text1 = text(0.45, 1.05, '�����ļ�ͷ�鿴��') ;
set(text1, 'FontSize', 13, 'Color', 'k' ,'FontWeight', 'bold', 'horizontalal', 'center');
% prepare to display header info on a GPR data file
top=0.9; left=0.05;
right=0.95; bottom=0.05;
label_height=0.07; spacing=0.02;
% Draw the text window frame
frame_border=0.02;
frame_position=[left-frame_border bottom-frame_border ...
    (right-left)+2*frame_border (top-bottom)+2*frame_border];
uicontrol( 'Style','frame', ...
    'Units','normalized', 'Position',frame_position, ...
    'BackgroundColor',[0.0 0.5 0.5]);
% Draw the text label
label_position=[left top-label_height (right-left) label_height];
uicontrol( 'Style','text', ...
    'Units','normalized', 'Position',label_position, ...
    'BackgroundColor',[0.0 0.5 0.5], 'ForegroundColor',[1 1 1], ...
    'String',title_string, 'fontsize', 10, 'fontweight', 'demi');
% Display the info box and display the info text
text_position=[left bottom (right-left) top-bottom-label_height-2*spacing];
uicontrol( 'Style','edit',...
    'tag', 'data_info_box', ...
    'Units','normalized', 'Max',18, 'String',data_info_text, ...
    'BackgroundColor',[1 1 1], ...
    'Position',text_position);

%% Print report
font_name = get(0, 'FixedWidthFontName');
set(findobj('tag','data_info_box'),...
    'String',str_data_info,...
    'horizontalal', 'left',...
    'fontname', font_name,...
    'fontsize', 18);
return
