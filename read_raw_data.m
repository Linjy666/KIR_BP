function [raw_data, data_info] = read_raw_data(file_name)
% 读取SFCW雷达数据，获取雷达参数信息
fid = fopen(file_name, 'rb');

data_info.ip_addr       = fread(fid, [1,16], '*char');% 雷达主机IP地址

% 波形参数
data_info.waveform      = fread(fid, 1, 'short');% 波形类型：0-步进频连续波SFCW；1-冲击脉冲
data_info.data_type     = fread(fid, 1, 'short');% 数据类型：0-原始数据ad；1-数字正交解调数据dpd；2-脉冲压缩数据rp
data_info.fl    = fread(fid, 1, 'float32');% SFCW信号的起始频率
data_info.fh    = fread(fid, 1, 'float32');% SFCW信号的终止频率
data_info.df    = fread(fid, 1, 'float32');% SFCW信号的频率步长
data_info.prf   = fread(fid, 1, 'float32');% 脉冲重复频率
data_info.fs	= fread(fid, 1, 'float32');% 接收机中频采样率
data_info.frequency_number	= fread(fid, 1, 'short');% SFCW频点数量
data_info.sample_number     = fread(fid, 1, 'short');% SFCW每个频点的中频采样样本数
data_info.ifft_number       = fread(fid, 1, 'short');% 脉冲压缩中逆傅里叶变换点数

% 阵列参数
data_info.channel_number    = fread(fid, 1, 'short');% 雷达系统全部通道数量=收发通道数量+参考通道数量
data_info.tr_pair_number    = fread(fid, 1, 'short');% 雷达系统收发通道数量
data_info.receiver_gain     = fread(fid, 1, 'short');% 雷达接收机增益
data_info.antenna_pos       = fread(fid, 400, 'float32');	% 收发天线组合坐标，按照tx,ty,tz,rx,ry,rz的顺序依次排列

data_info.reserved          = fread(fid, 6, 'short');       % 保留字节

data_info.antenna_pos = data_info.antenna_pos(1:data_info.tr_pair_number*3*2);
data_info.antenna_pos = reshape(data_info.antenna_pos, [data_info.tr_pair_number,3,2]);
if data_info.reserved(1) == 1
    data_info.frequency_number = data_info.frequency_number+1;
end

% SFCW
if data_info.waveform == 0
    switch (data_info.data_type)
        % 原始中频AD数据
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
                    msgbox(['第',num2str(abnormal_frame_index),'帧数据异常']);
                end
                data_info.frequency_number = data_info.frequency_number-1;
                sample_number_per_frame = data_info.sample_number*data_info.channel_number*data_info.frequency_number;
                raw_data = raw_data(1:sample_number_per_frame, :);
            end
            
        % 数字正交解调数据dpd
        case 1
            sample_number_per_frame = 2*data_info.channel_number*data_info.frequency_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'int32');
            frame_number = size(raw_data, 2);
            raw_data_iq = complex(raw_data(1:2:end,:), raw_data(2:2:end,:)); % 构建IQ两路复数据
            raw_data = reshape(raw_data_iq, [data_info.channel_number, data_info.frequency_number, frame_number]);
            raw_data = permute(raw_data, [2 1 3]);
            
            if data_info.reserved(1) == 1
                data_info.frequency_number = data_info.frequency_number-1;
                raw_data = raw_data(1:end-1,:,:);
            end
            
        % 通道校准数据rcc
        case 2
            sample_number_per_frame = 2*data_info.tr_pair_number*data_info.frequency_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'float32');
            if data_info.reserved(1) == 1
                data_info.frequency_number = data_info.frequency_number-1;
                raw_data = raw_data(1:data_info.tr_pair_number*data_info.frequency_number*2,:);
            end
            frame_number = size(raw_data, 2);
            raw_data_iq = complex(raw_data(1:2:end,:), raw_data(2:2:end,:)); % 构建IQ两路复数据
            raw_data = reshape(raw_data_iq, [data_info.tr_pair_number, data_info.frequency_number, frame_number]);
            raw_data = permute(raw_data, [2 1 3]);
            
        % 脉冲压缩数据rp
        case 3
            sample_number_per_frame = 2*data_info.tr_pair_number*data_info.ifft_number;
            raw_data = fread(fid, [sample_number_per_frame, inf], 'float32');
    end
end
fclose(fid);