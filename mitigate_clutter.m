function [result, backgorund] = mitigate_clutter(echo, backgorund, frame_count, mode, para)

switch mode
        % 脉冲对消
    case 'two_pulses_canceller'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            result = echo - backgorund;
            backgorund = echo;
        end
        % 脉冲对消
    case 'three_pulses_canceller'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = zeros([size(echo),2]);
            backgorund(:,:,1) = echo;
        elseif frame_count == 2
            result = zeros(size(echo));
            backgorund(:,:,2) = echo;
        else
            result = echo - 2*backgorund(:,:,2)+backgorund(:,:,1);
            backgorund(:,:,1) = backgorund(:,:,2);
            backgorund(:,:,2) = echo;
        end
        % 背景相消
    case 'background_subtraction'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            result = echo - backgorund;
            backgorund = echo/frame_count + backgorund*(1-1/frame_count);% 背景均值的估计
        end
        % 指数加权算法：alpha = (n-2)/(n-1)
    case 'exponential_weight'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            alpha = para;
            result = echo - backgorund;
            backgorund = alpha*backgorund+(1-alpha)*echo;
        end        
%         % 自适应背景相消
%     case 'adaptive_background_subtraction'
% %% 快速算法
%         % 分＜L和＞L两部分计算
%         % max(SCNR) = max((s-b)/b) = max(-b)
%         DstData(:, 1) = SrcData(:, 1); % 输出
%         C = ones(m, 1)*(1:Lmax); % 均值系数
%         Lopt(:, 1) = 0; % 最优抽头常数
%         B = zeros(m, Lmax); % 背景数据
%         B(:, 1) = SrcData(:, 1); % 背景数据初值
%         
%         for j = 2:Lmax
%             [T, I] = max(-B(:, 1:j-1), [], 2); % 计算出最优抽头长度，以及对应的背景数据
%             DstData(:, j) = SrcData(:, j) + T;
%             Lopt(:, j) = I;
%             B(:, 1:j) = [ SrcData(:, j),...
%                 ( B(:, 1:j-1).*C(:, 1:j-1) + SrcData(:, j)*ones(1, j-1) ) ./C(:, 2:j) ];% 更新背景数据
%         end
%         
%         for j = Lmax+1:n
%             [T, I] = max(-B(:, 3:end), [], 2);% 计算出最优抽头长度，以及对应的背景数据
%             DstData(:, j) = SrcData(:, j) + T;
%             Lopt(:, j) = I+2;
%             B = B - ( SrcData(:, j-(1:Lmax)) - SrcData(:, j)*ones(1, Lmax) ) ./C; % 更新背景数据
%         end

        % 滑动窗口背景相消
    case 'sliding_window'
        window_width = para;
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = zeros([size(echo),window_width]);
            backgorund(:,:,1) = echo;
        elseif frame_count <= window_width
            result = echo - mean(backgorund(:,:,1:frame_count-1),3);
            backgorund(:,:, frame_count) = echo;
        else
            result = echo - mean(backgorund(:,:,1:window_width), 3);
            backgorund(:,:, 1:window_width-1) = backgorund(:,:,2:window_width);
            backgorund(:,:, window_width) = echo;
        end
%         DstData(:, 1) = SrcData(:, 1); % 输出
%         B = SrcData(:, 1); % 背景数据初值
%         
%         for j = 2:Lmax
%             DstData(:, j) = SrcData(:, j) - B;
%             B = (B*(j-1) + SrcData(:, j)) / j;% 更新背景数据
%         end
%         
%         for j = Lmax+1:n
%             DstData(:, j) = SrcData(:, j) - B;
%             B = (B*Lmax - SrcData(:, j-Lmax) + SrcData(:, j) ) / Lmax; % 更新背景数据
%         end
end