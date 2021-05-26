function [x_tick, y_tick, z_tick, x_grid, y_grid, z_grid, index_mat, phase_mat] =...
    calculate_3d_image_grid(fl, dr, antenna_pos, image_info)
c = 3.0e8;
tr_pair_number = size(antenna_pos, 1);
tx = antenna_pos(:,1,1);ty = antenna_pos(:,2,1);tz = antenna_pos(:,3,1);
rx = antenna_pos(:,1,2);ry = antenna_pos(:,2,2);rz = antenna_pos(:,3,2);

x_min = image_info.x_min; dx = image_info.dx; x_max = image_info.x_max;
y_min = image_info.y_min; dy = image_info.dy; y_max = image_info.y_max;
z_min = image_info.z_min; dz = image_info.dz; z_max = image_info.z_max;
x_grid_size = (x_max-x_min)/dx+1; y_grid_size = (y_max-y_min)/dy+1; z_grid_size = (z_max-z_min)/dz+1;
x_tick = linspace(x_min, x_max, x_grid_size);
y_tick = linspace(y_min, y_max, y_grid_size);
z_tick = linspace(z_min, z_max, z_grid_size);
if strcmpi(image_info.axis_mode, 'polar')
    x_tick = x_tick*pi/180;
end
% 划分成像区域网格
[x_grid, y_grid, z_grid] = meshgrid(x_tick, y_tick, z_tick);
index_mat = zeros([size(x_grid), tr_pair_number]);
phase_mat = zeros([size(x_grid), tr_pair_number]);
% for i = 1:tr_pair_number
%     if strcmpi(image_info.axis_mode, 'xoy')
%         range = sqrt( (x_grid-rx(i)).^2 + (y_grid-ry(i)).^2 + rz(i).^2 ) +... % 接收天线与网格间距离
%             sqrt( (x_grid-tx(i)).^2 + (y_grid-ty(i)).^2 + tz(i).^2 );     % 发射天线与网格间距离
%     elseif strcmpi(image_info.axis_mode, 'polar')
%         range = sqrt( (y_grid.*sin(x_grid)-rx(i)).^2 + (y_grid.*cos(x_grid)-ry(i)).^2 + rz(i).^2 ) +... % 接收天线与网格间距离
%             sqrt( (y_grid.*sin(x_grid)-tx(i)).^2 + (y_grid.*cos(x_grid)-ty(i)).^2 + tz(i).^2 );     % 发射天线与网格间距离
%     end
%     index_mat(:,:,i) = round(range/(2*dr))+1;
%     phase_mat(:,:,i) = exp(1i*2*pi*fl*range/c); % 补偿的相位
% end
for i = 1:tr_pair_number
    if strcmpi(image_info.axis_mode, 'polar')
        r = sqrt(...
                    (y_grid.*cos(z_grid*pi/180).*sin(x_grid*pi/180)-tx(i)).^2 +...
                    (y_grid.*cos(z_grid*pi/180).*cos(x_grid*pi/180)-ty(i)).^2 +...
                    (y_grid.*sin(z_grid*pi/180)-tz(i)).^2 ...
                )...
                +...
            sqrt(...
                    (y_grid.*cos(z_grid*pi/180).*sin(x_grid*pi/180)-rx(i)).^2 +...
                    (y_grid.*cos(z_grid*pi/180).*cos(x_grid*pi/180)-ry(i)).^2 +...
                    (y_grid.*sin(z_grid*pi/180)-rz(i)).^2 ...
                );
    elseif strcmpi(image_info.axis_mode, 'xoy')
        r = sqrt(...
                    (x_grid-rx(i)).^2 + (y_grid-ry(i)).^2 + (z_grid-rz(i)).^2 ...
                )...
                +...
            sqrt(...
                    (x_grid-tx(i)).^2 + (y_grid-ty(i)).^2 + (z_grid-tz(i)).^2 ...
                );     % 发射天线与网格间距离
    end
    
%     if strcmpi(waveform_type, 'sfcw')
        phase_mat(:, :, :, i) = exp(1i*2*pi*fl*r/c); % 补偿的相位
%     end
    index_mat(:, :, :, i) = round(r/(2*dr))+1;
end

if strcmpi(image_info.axis_mode, 'polar')
    x_tick = x_tick*180/pi;
end


