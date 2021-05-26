function [bp_sum, bp_sum_cf] = back_projection_3d(range_profiles, index_mat, phase_mat)

%% BP成像
tr_pair_number = size(range_profiles, 2);
channel_images = zeros(size(index_mat));
for i = 1:tr_pair_number
    temp = range_profiles(:, i);
%     if strcmpi(waveform_type, 'sfcw')
        channel_images(:, :, :, i) = temp(index_mat(:,:,:,i)) .* phase_mat(:,:,:,i);
%     elseif strcmpi(waveform_type, 'impulse')
%         channel_images(:, :, :, i) = temp(index_mat(:,:,:,i));
%     end
end

% 原始BP
bp_sum = abs(sum(channel_images, 4));

% CF加权BP
cf = abs(bp_sum).^2 ./ (tr_pair_number*sum(abs(channel_images).^2, 4));
bp_sum_cf = bp_sum.*cf;