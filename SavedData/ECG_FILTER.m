clear all;
clc;
addpath(genpath('/home/ibagon/OpenBCI/OpenBCI_MATLAB/Matlab-Python/labstreaminglayer'))
% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an ECG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','ECG'); end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});
[vec,ts] = inlet.pull_sample();
start = ts;
eeg_record = [];
i = 1;
while ts - start < 1.5
    [vec,ts] = inlet.pull_sample();
    eeg_record = [eeg_record;vec];
    %fprintf('%d\t',i);
    %if i == 200
        %i=0;
    %end
    %fprintf('%.2f\t',vec);
    %fprintf('%.5f\n',ts);
    i = i+1;
    %vis_stream(vec);
end

%z=xlsread('xuewangnervous.xlsx','B7:B12000');
z=eeg_record[i]
[amp,time]=pan_tompkin(z,200,1);
hr=[];
time_new=[];
j=1;
for i=2:length(time)
    hr(j)=200*60/(time(i)-time(i-1));
    time_new(j)=time(i);
    j=j+1;
end
values=spcrv([[time_new(1) time_new time_new(end)];[hr(1) hr hr(end)]],3);
figure;
plot(values(1,:),values(2,:));
title('Heart Rate Calculated');
%{
rate=z(100);
sig=z;
lensig=length(sig);
wtsig1=cwt(sig,6,'mexh');
lenwtsig1=length(wtsig1);
%wtsig1(1:8)=0;
%wtsig1(lenwtsig1-10:lenwtsig1)=0;
y=wtsig1;
yabs=abs(y);
sigtemp=y;
siglen=length(y);
sigmax=[];
for i=1:siglen-2
    if (y(i+1)>y(i)&&y(i+1)>y(i+2))||(y(i+1)<y(i)&&y(i+1)<y(i+2))
        sigmax=[sigmax;sigtemp(i+1),i+1];
    end
end

%取阈值,阈值为相对幅值的差的60%
%set threshold 
thrtemp=sort(sigmax);
thrlen=length(sigmax);
thr=0;
for i=(thrlen-1000):thrlen
    thr=thr+thrtemp(i);
end
thrmax=thr/1000;               %最大幅度平均值，1000个最大幅值点的平均值 average of max_amptitude
zerotemp=sort(y);
zerovalue=0;
for i=1:1000
    zerovalue=zerovalue+zerotemp(i);
end
zerovalue=zerovalue/1000;    %最小幅度平均值，对消幅度，1000个最小幅值点的平均值 average of min_amptitude
 
thr=(thrmax-zerovalue)*0.70; %最大、最小幅度的差值的70%为判别R波的阈值  thresh=0.7*(max-min)                    
 
%定位R波
%find location of R wave
rvalue=[];
for i=1:thrlen
    if sigmax(i,1)>thr
        rvalue=[rvalue;sigmax(i,2)];
    end
end
rvalue_1=rvalue;

%计算心率
%calculate heart rate
hr=[];
time_new=[];
j=1;
for i=2:length(rvalue_1)
    if (rvalue_1(i)-rvalue(i-1))>120
        hr(j)=200*60/(rvalue_1(i)-rvalue(i-1));
        time_new(j)=rvalue_1(i);
        j = j +1;
    end
end

%排除误检，如果相邻两个极大值间距小于0.4?，则去掉幅度较小的一个
%eliminate 'fake' R wave
lenvalue=length(rvalue);
i=2;
while i<=lenvalue
      if (rvalue(i)-rvalue(i-1))*rate<0.4
          if yabs(rvalue(i))>yabs(rvalue(i-1))
              rvalue(i-1)=[];
          else
              rvalue(i)=[];
          end
          lenvalue=length(rvalue);
          i=i-1;
      end
      i=i+1;
end        
 
lenvalue=length(rvalue);

%打印
%draw plots
figure(1);
subplot(3,1,1),plot(sig);
title('ECG Raw Data');
subplot(3,1,2),plot(1:lensig,wtsig1,rvalue_1,wtsig1(rvalue_1),'r.');
title('ECG Filtered Data with R-Wave Pattern Location Stamped');
values=spcrv([[time_new(1) time_new time_new(end)];[hr(1) hr hr(end)]],3);
subplot(3,1,3),plot(values(1,:),values(2,:));
title('Heart Rate Calculated');
%}

