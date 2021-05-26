function [p, y, m1, n1, m2, n2, min_index, max_index] = evaluate_linear_degree(channel_data)

p = phase(channel_data);
% �����λ�ֲ�����
x = (1:length(p))';
A = polyfit(x,p,1);
y = polyval(A,x);
% �������Զ�
d1 = (p-y)./(max(y)-min(y));
[m1, min_index] = min(d1); [n1, max_index] = max(d1);

d2 = p-y;
m2 = min(d2); n2 = max(d2);

% figure,
% plot(p,'-');
% grid on;title('��λ�ֲ�');
% hold on, plot(y, 'r-');
% xlabel('Ƶ��');ylabel('��λ');xlim([0 frequency_number]);
% legend('ԭʼ��λ�ֲ�', '�����λ�ֲ�');
% annotation(gcf,'textbox',...
%     [0.444 0.64 0.44 0.138],...
%     'Color',[1 0 0],...
%     'String',{['[',num2str(100*m, '%6.2f'),'%, ', num2str(100*n, '%6.2f'),'%]']},...
%     'LineStyle','none',...
%     'FontWeight','bold',...
%     'FontSize',24,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off');
end