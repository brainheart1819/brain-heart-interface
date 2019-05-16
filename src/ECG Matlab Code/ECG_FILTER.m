% % clear all;
% % clc;
% % z=xlsread('kuicao_tender.xlsx','B07:B180000');
% % ramp=zeros();
% % ramp_raw=zeros();
% % samp=zeros();
% % samp_raw=zeros();
% % qamp=zeros();
% % qamp_raw=zeros();
% % tamp=zeros();
% % tamp_raw=zeros();
% % [rlocation,amp,ecg_h]=pan_tompkin(z,200,1);
% % for i=1:length(rlocation)
% %     ramp(i) = ecg_h(rlocation(i));
% %     ramp_raw(i) = z(rlocation(i));
% % end
% % [slocation]=SWavedetection(rlocation,ecg_h);
% % for i=1:length(slocation)
% %     samp(i) = ecg_h(slocation(i));
% %     samp_raw(i) = z(slocation(i));
% % end
% % [qlocation]=QWavedetection(rlocation,ecg_h);
% % for i=1:length(qlocation)
% %     qamp(i) = ecg_h(qlocation(i));
% %     qamp_raw(i) = z(qlocation(i));
% % end
% % [tlocation]=TWavedetection(slocation,qlocation,ecg_h);
% % for i=1:length(tlocation)
% %     tamp(i) = ecg_h(tlocation(i));
% %     tamp_raw = z(tlocation(i));
% % end
% % 
% % plot(z/200);
% % title('kuicaotender');
% % axis ([0 2000 -5 5]);
% % hold on,plot(ecg_h);
% % hold on,scatter(rlocation,ramp,'b','filled');
% % hold on,scatter(slocation,samp,'r','filled');
% % hold on,scatter(qlocation,qamp,'g','filled');
% % hold on,scatter(tlocation,tamp,'m','filled');
% % legend('raw data','filtered signal','R wave peak location',...
% %     'S wave peak location','Q wave peak location','T wave peak location');
% % 
% % %1.mean of all R-wave in one minute
% % mean_r=mean(ramp_raw);
% % 
% % %2.R_range of all R_wave in one minute
% % range_r=max(ramp_raw)-min(ramp_raw);
% % 
% % %3.the mean of all Q-wave
% % mean_q=mean(qamp_raw);
% % 
% % %4.the standard deviation of all S-wave
% % std_s=std(samp_raw);
% % 
% % %5.the mean of each S-T interval
% % st_interval=zeros();
% % if (length(slocation) > length(tlocation))
% %     for i=1:length(tlocation)
% %         st_interval(i) = tlocation(i) - slocation(i);
% %     end
% % elseif (length(slocation) < length(tlocation))
% %     for i=1:length(slocation)
% %         st_interval(i) = tlocation(i+1) - slocation(i);
% %     end
% % else
% %     if (slocation(1) > tlocation(1))
% %         for i=1:(length(slocation)-1)
% %             st_interval(i) = tlocation(i+1) - slocation(i);
% %         end
% %     else
% %         for i=1:length(slocation)
% %             st_interval(i) = tlocation(i) - slocation(i);
% %         end
% %     end
% % end
% % mean_st = mean(st_interval);
% % 
% % %6.the mean of heart rate calculated (take advantage of RR interval)
% % hr=zeros();
% % for i=2:length(rlocation)
% %     hr(i)=200*60/(rlocation(i)-rlocation(i-1));
% % end
% % mean_hr=mean(hr);
% % 
% % %7.the standard deviation of all heart rates
% % std_hr=std(hr);
% % 
% % %8.the mean of all heart rate variability calculated (a common way in which
% % %to analyse HRV is a time-domain method called RMSSD. This is the Root Mean
% % %Square of Successive Differenves between wach heartbeat)
% % rrinterval_square=zeros();
% % for i=2:length(rlocation)
% %     rrinterval_square(i)=(rlocation(i)-rlocation(i-1))^2;
% % end
% % rmssd=mean(rrinterval_square(2:end));
% % 
% % %9.the mean of RR interval
% % rrinterval=zeros();
% % for i=2:length(rlocation)
% %     rrinterval(i)=rlocation(i)-rlocation(i-1);
% % end
% % mean_rr_interval=mean(rrinterval(2:end));
% % 
% % %10.the mean of all the ECG raw signal data
% % mean_raw=mean(z);
% % 
% % 
% % 
% % %frequency domain analyse
% % [Pxx,F] = periodogram(z,[],numel(z),1);
% % %1.the mean of the power spectral density of very low frequency data
% % vlf = bandpower(Pxx,F,[0 0.04],'psd');
% % 
% % %2.the mean of the power spectral density of low frequency data
% % lf = bandpower(Pxx,F,[0.04 0.15],'psd');
% % 
% % %3.the mean of the power spectral density of high frequency data
% % hf = bandpower(Pxx,F,[0.15 0.4],'psd');
% % 
% % totalpower = bandpower(Pxx,F,'psd');
% % 
% % %4.the percentage power of low frquency data occupied in total power
% % %spectrum expect vlf
% % lf_ratio = lf/(totalpower-vlf);
% % 
% % %5.the percentage power of high frquency data occupied in total power
% % %spectrum expect vlf
% % hf_ratio = hf/(totalpower-vlf);
% % 
% % %6.the ratio of low frquency data psd versus high frequency data psd
% % lf_hf_ratio = lf/hf;




clear all;
clc;
z=xlsread('0324data/xuezhiwang_tender.xlsx','B:B');
ramp=zeros();
ramp_raw=zeros();
samp=zeros();
samp_raw=zeros();
qamp=zeros();
qamp_raw=zeros();
tamp=zeros();
tamp_raw=zeros();
[rlocation,amp,ecg_h]=pan_tompkin(z,200,1);
for i=1:length(rlocation)
    ramp(i) = ecg_h(rlocation(i));
    ramp_raw(i) = z(rlocation(i));
end
[slocation]=SWavedetection(rlocation,ecg_h);
for i=1:length(slocation)
    samp(i) = ecg_h(slocation(i));
    samp_raw(i) = z(slocation(i));
end
[qlocation]=QWavedetection(rlocation,ecg_h);
for i=1:length(qlocation)
    qamp(i) = ecg_h(qlocation(i));
    qamp_raw(i) = z(qlocation(i));
end
[tlocation]=TWavedetection(slocation,qlocation,ecg_h);
for i=1:length(tlocation)
    tamp(i) = ecg_h(tlocation(i));
    tamp_raw = z(tlocation(i));
end
plot(z/200);
axis ([0 18000 -5 5]);
hold on,plot(ecg_h);
hold on,scatter(rlocation,ramp,'b','filled');
hold on,scatter(slocation,samp,'r','filled');
hold on,scatter(qlocation,qamp,'g','filled');
hold on,scatter(tlocation,tamp,'m','filled');
legend('raw data','filtered signal','R wave peak location',...
    'S wave peak location','Q wave peak location','T wave peak location');

%1.mean of all R-wave in one minute
mean_r=mean(ramp_raw);

mean_r_buffer=zeros();
mean_r_array=zeros();
j = 1;
k = 1;
l = 1;
for i=1:length(z)
    if (i > j*1000)
        mean_r_array(j) = mean(mean_r_buffer);
        mean_r_buffer = zeros();
        j = j + 1;
        k = 1;
    elseif (i == length(z) || l >= length(rlocation))
        mean_r_array(j) = mean(mean_r_buffer);
        break;
    elseif (i == rlocation(l))
        mean_r_buffer(k) = z(rlocation(l));
        k = k + 1;
        l = l + 1;
    end
    
end



%2.R_range of all R_wave in one minute
range_r=max(ramp_raw)-min(ramp_raw);

%3.the mean of all Q-wave
mean_q=mean(qamp_raw);

mean_q_buffer=zeros();
mean_q_array=zeros();
j = 1;
k = 1;
l = 1;
for i=1:length(z)
    if (i > j*1000)
        mean_q_array(j) = mean(mean_q_buffer);
        mean_q_buffer = zeros();
        j = j + 1;
        k = 1;
    elseif (i == length(z) || l >= length(qlocation))
        mean_q_array(j) = mean(mean_q_buffer);
        break;
    elseif (i == qlocation(l))
        mean_q_buffer(k) = z(qlocation(l));
        k = k + 1;
        l = l + 1;
    end
end

%** mean_r/mean_q
r_q_ratio = zeros();
for i = 1:min(length(mean_r_array), length(mean_q_array))
    r_q_ratio(i) = abs(mean_r_array(i) - mean_r) / abs(mean_q_array(i) - mean_q);
end

%4.the standard deviation of all S-wave
std_s=std(samp_raw);

%5.the mean of each S-T interval
st_interval=zeros();
if (length(slocation) > length(tlocation))
    for i=1:length(tlocation)
        st_interval(i) = tlocation(i) - slocation(i);
    end
elseif (length(slocation) < length(tlocation))
    for i=1:length(slocation)
        st_interval(i) = tlocation(i+1) - slocation(i);
    end
else
    if (slocation(1) > tlocation(1))
        for i=1:(length(slocation)-1)
            st_interval(i) = tlocation(i+1) - slocation(i);
        end
    else
        for i=1:length(slocation)
            st_interval(i) = tlocation(i) - slocation(i);
        end
    end
end
mean_st = mean(st_interval);

%6.the mean of heart rate calculated (take advantage of RR interval)
hr=zeros();
for i=2:length(rlocation)
    hr(i)=200*60/(rlocation(i)-rlocation(i-1));
end
mean_hr=mean(hr);

%7.the standard deviation of all heart rates
std_hr=std(hr);

std_hr_buffer=zeros();
std_hr_array=zeros();
j = 1;
k = 1;
l = 1;
for i=1:length(z)
    if (i > j*1000)
        std_hr_array(j) = std(std_hr_buffer);
        std_hr_buffer = zeros();
        j = j + 1;
        k = 1;
    elseif (i == length(hr) || l >= length(rlocation))
        std_hr_array(j) = std(std_hr_buffer);
    elseif (i == rlocation(l))
        std_hr_buffer(k) = hr(l);
        k = k + 1;
        l = l + 1;
    end
end




%8.the mean of all heart rate variability calculated (a common way in which
%to analyse HRV is a time-domain method called RMSSD. This is the Root Mean
%Square of Successive Differenves between wach heartbeat)
rrinterval_square=zeros();
for i=2:length(rlocation)
    rrinterval_square(i)=(rlocation(i)-rlocation(i-1))^2;
end
rmssd=mean(rrinterval_square(2:end));

%9.the mean of RR interval
rrinterval=zeros();
for i=2:length(rlocation)
    rrinterval(i)=rlocation(i)-rlocation(i-1);
end
mean_rr_interval=mean(rrinterval(2:end));

%10.the mean of all the ECG raw signal data
mean_raw=mean(z);

mean_raw_buffer=zeros();
mean_raw_array=zeros();
j = 1;
k = 1;
for i=1:length(z)
    if (i > j*1000)
        mean_raw_array(j) = mean(mean_raw_buffer);
        mean_raw_buffer = zeros();
        j = j + 1;
        k = 1;
    elseif (i == length(z))
        mean_raw_array(j) = mean(mean_raw_buffer);
    end
    mean_raw_buffer(k) = z(i);
    k = k + 1;
end

n = min(length(r_q_ratio), min(length(std_hr_array), length(mean_raw_array)));
r_q_ratio = r_q_ratio(1:n);
std_hr_array = std_hr_array(1:n);
mean_raw_array = mean_raw_array(1:n);
output = [r_q_ratio;std_hr_array;mean_raw_array];
csvwrite("0324data/xuezhiwangtender.csv",output.')



%frequency domain analyse
[Pxx,F] = periodogram(z,[],numel(z),1);
%1.the mean of the power spectral density of very low frequency data
vlf = bandpower(Pxx,F,[0 0.04],'psd');

%2.the mean of the power spectral density of low frequency data
lf = bandpower(Pxx,F,[0.04 0.15],'psd');

%3.the mean of the power spectral density of high frequency data
hf = bandpower(Pxx,F,[0.15 0.4],'psd');

totalpower = bandpower(Pxx,F,'psd');

%4.the percentage power of low frquency data occupied in total power
%spectrum expect vlf
lf_ratio = lf/(totalpower-vlf);

%5.the percentage power of high frquency data occupied in total power
%spectrum expect vlf
hf_ratio = hf/(totalpower-vlf);

%6.the ratio of low frquency data psd versus high frequency data psd
lf_hf_ratio = lf/hf;


fprintf('DONE!\n');




%{
addpath(genpath('/home/ibagon/OpenBCI/OpenBCI_MATLAB/Matlab-Python/labstreaminglayer'))

% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an ECG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG');
end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

ecg = zeros();
[vec,ts] = inlet.pull_sample();
start = ts;
while ts-start<15
    [vec,ts] = inlet.pull_sample();
    fprintf('%.2f\n',vec(1));
end

while ts-start<=4
    %[amp,time,ecg_h] = pan_tompkin(ecg,200,1);
    %plot(ecg);
end
    %i = i+1;
    %j = j+0.5;
    %addpoints(h1,i,vec(1));
    %drawnow
    %title('raw data');
    %axis([j j+200 -500 150]);
    %count=0;
    %while rem(i,4000)==0
        %[amp,time,ecg_h]=pan_tompkin(ecg,200,1); 
        %hr=[];
        %time_new=[];
        %k=1;
        %for i=2:length(time)
            %hr(i)=200*60/(time(i)-time(i-1));
            %count = count+1;
            %time_new(k)=time(i);
            %k=k+1;
        %end
            %values=spcrv([[time_new(1) time_new time_new(end)];[hr(1) hr hr(end)]],3);
            %figure;
            %plot(values(1,:),values(2,:));
            %title('Heart Rate Calculated');
            %figure,az(1)=plot(ecg_h);title('QRS on Filtered Signal');axis tight;
            %hold on,scatter(time,amp,'m');
    %end
%end
        %{
        values=spcrv([[time_new(1) time_new time_new(end)];[hr(1) hr hr(end)]],3);
        figure;
        t=[0];
        beat=[0];
        p=plot(t,beat,'EraseMode','background','MarkerSize',5);
        grid on;
        %for i=1:1000
            t=values(1,:);
            beat=values(2,:);
            set(p,'XData',t,'YData',beat);
            drawnow;
            pause(0.1);
        %end
        drawnow;
        title('Heart Rate Calculated');
       %fprintf('%.2f\n',amp);
       i=1;
       ecg = [];
        %}
%}

  %{  
    hr=zeros();
    for i=2:length(time)-1
    hr=200*60/(time(i)-time(i-1));
    end
    figure;
    scatter(time,hr);
    title('Heart Rate Calculated');
%}    
