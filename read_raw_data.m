function [raw_data, data_info] = read_raw_data(file_name)
% ��ȡSFCW�״����ݣ���ȡ�״������Ϣ
fid = fopen(file_name, 'rb');

data_info.ip_addr       = fread(fid, [1,16], '*char');% �״�����IP��ַ

% ���β���
data_info.waveform      = fread(fid, 1, 'short');% �������ͣ�0-����Ƶ������SFCW��1-�������
data_info.data_type     = fread(fid, 1, 'short');% �������ͣ�0-ԭʼ����ad��1-���������������dpd��2-����ѹ������rp
data_info.fl    = fread(fid, 1, 'float32');% SFCW�źŵ���ʼƵ��
data_info.fh    = fread(fid, 1, 'float32');% SFCW�źŵ���ֹƵ��
data_info.df    = fread(fid, 1, 'float32');% SFCW�źŵ�Ƶ�ʲ���
data_info.prf   = fread(fid, 1, 'float32');% �����ظ�Ƶ��
data_info.fs	= fread(fid, 1, 'float32');% ���ջ���Ƶ������
data_info.frequency_number	= fread(fid, 1, 'short');% SFCWƵ������
data_info.sample_number     = fread(fid, 1, 'short');% SFCWÿ��Ƶ�����Ƶ����������
data_info.ifft_number       = fread(fid, 1, 'short');% ����ѹ�����渵��Ҷ�任����

% ���в���
data_info.channel_number    = fread(fid, 1, 'short');% �״�ϵͳȫ��ͨ������=�շ�ͨ������+�ο�ͨ������
data_info.tr_pair_number    = fread(fid, 1, 'short');% �״�ϵͳ�շ�ͨ������
data_info.receiver_gain     = fread(fid, 1, 'short');% �״���ջ�����
data_info.antenna_pos       = fread(fid, 400, 'float32');	% �շ�����������꣬����tx,ty,tz,rx,ry,rz��˳����������

data_info.reserved          = fread(fid, 6, 'short');       % �����ֽ�

data_info.antenna_pos = data_info.antenna_pos(1:data_info.tr_pair_number*3*2);
data_info.antenna_pos = reshape(data_info.antenna_pos, [data_info.tr_pair_number,3,2]);
if data_info.reserved(1) == 1
    data_info.frequency_number = data_info.frequency_number+1;
end

% SFCW
if data_info.waveform == 0
    switch (data_info.data_type)
        % ԭʼ��ƵAD����
        case 0
            sample_number_per_frame = data_info.sample_number*data_info.channel_number*data_info.frequency_number;
%             raw_data = fread(fid, 2, 'int16');
            raw_data = fread(fid, [sample_number_per_frame, inf], 'int16');
%             raw_data = raw_data(:,1:end-1);
            
            if data_info.reserved(1) == 1
                frame_number = size(raw_data, 2);
                right_flag = mod(0:frame_number-1, 2^16-1)-1;
                data_flag = raw_data(sample_number_per_frame-data_info.sample_number*data_info.channel_number+1:end, :);
                compare_error = mean(data_flag(1:data_info.sample_number, :))+...
                    mean(data_flag(data_info.sample_number+1:end, :))+...
                    -right_flag;
                abnormal_frame_index = find(compare_error~=0, 1);
                if ~isempty(abnormal_frame_index)
                    msgbox(['��',num2str(abnormal_frame_index),'֡�����쳣']);
                end
                data_info.frequency_number = data_info.frequency_number-1;
                sample_number_per_frame = data_info.sample_number*data_info.channel_number*data_info.frequency_number;
                raw_data = raw_data(1:sample_number_per_frame, :);
            end
            
        % ���������������dpd
        case 1
            sample_number_per_frame = 2*data_info.channel_number*data_info.frequency_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'int32');
            frame_number = size(raw_data, 2);
            raw_data_iq = complex(raw_data(1:2:end,:), raw_data(2:2:end,:)); % ����IQ��·������
            raw_data = reshape(raw_data_iq, [data_info.channel_number, data_info.frequency_number, frame_number]);
            raw_data = permute(raw_data, [2 1 3]);
            
            if data_info.reserved(1) == 1
                data_info.frequency_number = data_info.frequency_number-1;
                raw_data = raw_data(1:end-1,:,:);
            end
            
        % ͨ��У׼����rcc
        case 2
            sample_number_per_frame = 2*data_info.tr_pair_number*data_info.frequency_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'float32');
            if data_info.reserved(1) == 1
                data_info.frequency_number = data_info.frequency_number-1;
                raw_data = raw_data(1:data_info.tr_pair_number*data_info.frequency_number*2,:);
            end
            frame_number = size(raw_data, 2);
            raw_data_iq = complex(raw_data(1:2:end,:), raw_data(2:2:end,:)); % ����IQ��·������
            raw_data = reshape(raw_data_iq, [data_info.tr_pair_number, data_info.frequency_number, frame_number]);
            raw_data = permute(raw_data, [2 1 3]);
            
        % ����ѹ������rp
        case 3
            sample_number_per_frame = 2*data_info.tr_pair_number*data_info.ifft_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'float32');
    end
end
fclose(fid);