ch1=squeeze(handles.radar_data(:,4,:));
load('bp_filter_0.1_0.15_4_4.05_15.1994.mat');
ch1_flt = filter(bp_filter,1,ch1,[],2);
n = 512;
s1 = abs(fft(ch1_flt,n,2));

fs = 15.1994;
f_tick=(0:n/2-1)/n*fs;
ru = 150;
ifft_number = 4096;
dr = ru/ifft_number;
y_min = 5;
y_max = 30;
y_min_index = round(y_min/dr)+1;
y_max_index = round(y_max/dr)+1;
index = y_min_index:y_max_index;

y_tick=index*dr;
figure,imagesc(f_tick, y_tick, abs(s1(index, 1:n/2)));grid on
colormap('jet');set(gca, 'YDir', 'normal');
xlabel('ÆµÂÊ(Hz)');ylabel('¾àÀë(m)');