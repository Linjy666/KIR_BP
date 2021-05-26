function [result, backgorund] = mitigate_clutter(echo, backgorund, frame_count, mode, para)

switch mode
        % �������
    case 'two_pulses_canceller'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            result = echo - backgorund;
            backgorund = echo;
        end
        % �������
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
        % ��������
    case 'background_subtraction'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            result = echo - backgorund;
            backgorund = echo/frame_count + backgorund*(1-1/frame_count);% ������ֵ�Ĺ���
        end
        % ָ����Ȩ�㷨��alpha = (n-2)/(n-1)
    case 'exponential_weight'
        if frame_count == 1
            result = zeros(size(echo));
            backgorund = echo;
        else
            alpha = para;
            result = echo - backgorund;
            backgorund = alpha*backgorund+(1-alpha)*echo;
        end        
%         % ����Ӧ��������
%     case 'adaptive_background_subtraction'
% %% �����㷨
%         % �֣�L�ͣ�L�����ּ���
%         % max(SCNR) = max((s-b)/b) = max(-b)
%         DstData(:, 1) = SrcData(:, 1); % ���
%         C = ones(m, 1)*(1:Lmax); % ��ֵϵ��
%         Lopt(:, 1) = 0; % ���ų�ͷ����
%         B = zeros(m, Lmax); % ��������
%         B(:, 1) = SrcData(:, 1); % �������ݳ�ֵ
%         
%         for j = 2:Lmax
%             [T, I] = max(-B(:, 1:j-1), [], 2); % ��������ų�ͷ���ȣ��Լ���Ӧ�ı�������
%             DstData(:, j) = SrcData(:, j) + T;
%             Lopt(:, j) = I;
%             B(:, 1:j) = [ SrcData(:, j),...
%                 ( B(:, 1:j-1).*C(:, 1:j-1) + SrcData(:, j)*ones(1, j-1) ) ./C(:, 2:j) ];% ���±�������
%         end
%         
%         for j = Lmax+1:n
%             [T, I] = max(-B(:, 3:end), [], 2);% ��������ų�ͷ���ȣ��Լ���Ӧ�ı�������
%             DstData(:, j) = SrcData(:, j) + T;
%             Lopt(:, j) = I+2;
%             B = B - ( SrcData(:, j-(1:Lmax)) - SrcData(:, j)*ones(1, Lmax) ) ./C; % ���±�������
%         end

        % �������ڱ�������
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
%         DstData(:, 1) = SrcData(:, 1); % ���
%         B = SrcData(:, 1); % �������ݳ�ֵ
%         
%         for j = 2:Lmax
%             DstData(:, j) = SrcData(:, j) - B;
%             B = (B*(j-1) + SrcData(:, j)) / j;% ���±�������
%         end
%         
%         for j = Lmax+1:n
%             DstData(:, j) = SrcData(:, j) - B;
%             B = (B*Lmax - SrcData(:, j-Lmax) + SrcData(:, j) ) / Lmax; % ���±�������
%         end
end