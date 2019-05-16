% function my_eeg_signal_demo()
%% Initialize the variables
dirPath = "F:\data\3_24_set";
% dirPath = "F:\data\10_25_set";
% dirPath = "F:\Brain_Heart_Computer_Interface\313_set";
nameSet = [ "Kuicao" "Xuezhi"];
% nameSet = ["Kuicao","Pengwu","Xuezhiwang","Xuewang"];
% nameSet = ["Kuicao","Xuewang"];
% emotionSet = ["Happy" "Tender" "Tension" "Sad"];
emotionSet = ["Happy" "Sad" "Tender"];
data = [];
fs = 200;
channel_num = 2;
figure_num = 1;
%% read all EEG signal data from set files
for name = nameSet
    for emotion = emotionSet
        fileName = strcat(name, "_", emotion,".set");
        EEG = pop_loadset(char(fileName),char(dirPath));
        tempData = EEG.data(3:4,60*fs:150*fs);
        data = [data; tempData];
    end
end
        
% %% Calculate the statistic feature for raw dataset
M = mean(data,2);
S = std(data,0,2);
M1 = mean(abs(diff(data,1,2)),2);
S1 = M1./S;
M2 = mean(abs(diff(data,2,2)),2);
S2 = M2./S;
% 
% for i = 1:length(nameSet)
%     for j = 1:length(emotionSet) * channel_num
%         M_final(i,j) = M((i-1) * length(emotionSet) * channel_num + j);
%         S_final(i,j) = S((i-1) * length(emotionSet) * channel_num + j);
%         M1_final(i,j) = M1((i-1) * length(emotionSet) * channel_num + j);
%         S1_final(i,j) = S1((i-1) * length(emotionSet) * channel_num + j);
%         M2_final(i,j) = M2((i-1) * length(emotionSet) * channel_num + j);
%         S2_final(i,j) = S2((i-1) * length(emotionSet) * channel_num + j);
%     end
% end
% 
% figure(figure_num);
% figure_num = figure_num + 1;
% subplot(3,2,1);
% my_statistic_plot(M_final,8,'Mean Value',nameSet,emotionSet);
% subplot(3,2,2);
% my_statistic_plot(S_final,8,'Standard Variance',nameSet,emotionSet);
% subplot(3,2,3);
% my_statistic_plot(M1_final,8,'First Difference',nameSet,emotionSet);
% subplot(3,2,4);
% my_statistic_plot(S1_final,8,'Normalized First Difference',nameSet,emotionSet);
% subplot(3,2,5);
% my_statistic_plot(M2_final,8,'Second Difference',nameSet,emotionSet);
% subplot(3,2,6);
% my_statistic_plot(S2_final,8,'Normalized Second Difference',nameSet,emotionSet);
% 
% csvwrite('mean_value.csv', M_final);
% csvwrite('standard_variance.csv', S_final);
% csvwrite('first difference.csv', M1_final);
% csvwrite('normalized_first_difference.csv', S1_final);
% csvwrite('second_difference.csv', M2_final);
% csvwrite('normalized_second_difference.csv', S2_final);
      
%% denoise the data
data_without_BL = data - M;
filtSpec.order = 10;
filtSpec.range = [1,50];
filtPts = fir1(filtSpec.order, 2/fs*filtSpec.range);
filteredData = filter(filtPts, 1, data_without_BL, [], 2);
for i = 1:size(filteredData,1)
    filteredData(i,:) = reconstruct_data(filteredData(i,:));
end

%% Calculate the PSD for different bandwidth
BP = [];
BP_final = [];
BW = ["delta" "theta" "alpha" "beta" "low_beta" "high_beta" "gamma"];
for i = 1:size(filteredData,1)
    BP(:,:,i) = cal_BP(filteredData(i,:),fs);
end

for i = 1:size(BP,1)
    BP_final(:,:,i) = squeeze(BP(i,:,:))';
    csvwrite(char(strcat(BW(i),"_BP.csv")),BP_final(:,:,i))
end
BP_final2 = [];
for i = 1:size(BP_final,1)/2
    BP_final2(i,:,:) = squeeze(mean(BP_final(2*i-1:2*i,:,:)));
end
% BP_BASE1 = mean(squeeze(BP_final2(3,:,:)),1);
% BP_BASE2 = mean(squeeze(BP_final2(6,:,:)),1);
% for i=1:2
%     for j=1:139
%         BP_final2(i,j,:) = squeeze(BP_final2(i,j,:)) - BP_BASE1';
%     end
% end
% 
% for i=4:5
%     for j=1:139
%         BP_final2(i,j,:) = squeeze(BP_final2(i,j,:)) - BP_BASE2';
%     end
% end

for bw = 1:length(BW)    
    figure(figure_num)
    figure_num = figure_num + 1;
    for i = 1:length(nameSet)
        subplot(2,2,i);
        p = plot(squeeze(BP_final(8*(i-1)+1:8*i,:,bw))');
        p(1).Color = 'red';
        p(2).Color = 'red';
        p(3).Color = 'blue';
        p(4).Color = 'blue';
        p(5).Color = 'green';
        p(6).Color = 'green';
        p(7).Color = 'yellow';
        p(8).Color = 'yellow';
        legend([p(1) p(3) p(5) p(7)], emotionSet);
        title(char(strcat(nameSet(i)," ", BW(bw)," BP")));
    end
end
    
    
Alpha_beta_ratio =  BP_final(:,:,3)./BP_final(:,:,4);
Low_beta_high_beta_ratio = BP_final(:,:,5)./BP_final(:,:,6);

figure(figure_num)
figure_num = figure_num + 1;
for i = 1:length(nameSet)
    subplot(2,2,i);
    p = plot(squeeze(Alpha_beta_ratio(8*(i-1)+1:8*i,:))');
    p(1).Color = 'red';
    p(2).Color = 'red';
    p(3).Color = 'blue';
    p(4).Color = 'blue';
    p(5).Color = 'green';
    p(6).Color = 'green';
    p(7).Color = 'yellow';
    p(8).Color = 'yellow';
    legend([p(1) p(3) p(5) p(7)], emotionSet);
    title(char(strcat(nameSet(i)," ratio of alpha and beta signal")));
end

figure(figure_num)
figure_num = figure_num + 1;
for i = 1:length(nameSet)
    subplot(2,2,i);
    p = plot(squeeze(Low_beta_high_beta_ratio(8*(i-1)+1:8*i,:))');
    p(1).Color = 'red';
    p(2).Color = 'red';
    p(3).Color = 'blue';
    p(4).Color = 'blue';
    p(5).Color = 'green';
    p(6).Color = 'green';
    p(7).Color = 'yellow';
    p(8).Color = 'yellow';
    legend([p(1) p(3) p(5) p(7)], emotionSet);
    title(char(strcat(nameSet(i)," ratio of low beta and high beta signal")));
end

BP_Mean = [];
for i = 1:size(BP_final,3)
    BP_Mean(:,i) = mean(squeeze(BP_final(:,:,i)),2); 
end

for bw = 1:length(BW)    
    figure(figure_num)
    figure_num = figure_num + 1;
    for i = 1:length(nameSet)
        subplot(2,2,i);
        p = plot((squeeze(BP_final(8*(i-1)+1:8*i,:,bw))./BP_Mean(8*(i-1)+1:8*i,bw))');
        p(1).Color = 'red';
        p(2).Color = 'red';
        p(3).Color = 'blue';
        p(4).Color = 'blue';
        p(5).Color = 'green';
        p(6).Color = 'green';
        p(7).Color = 'yellow';
        p(8).Color = 'yellow';
        legend([p(1) p(3) p(5) p(7)], emotionSet);
        title(char(strcat(nameSet(i)," ", BW(bw)," BP ratio")));
    end
end
%% Alpha Asymmetry
alpha_asys = [];
for i = 1:size(squeeze(PSD_final(:,:,3)),1)/2
    alpha_asys = [alpha_asys; squeeze(PSD_final(2*i-1,:,3)./PSD_final(2*i,:,3))];
end
alpha_asys = log(alpha_asys);


figure(figure_num)
figure_num = figure_num + 1;
for i = 1:length(nameSet)
    subplot(2,2,i);
    p = plot(squeeze(alpha_asys(4*(i-1)+1:4*i,:))');
    p(1).Color = 'red';
    p(2).Color = 'blue';
    p(3).Color = 'green';
    p(4).Color = 'yellow';
    legend(emotionSet);
    title(char(strcat(nameSet(i)," alpha asymmetry value")));
end



%% PLV feature
epoch = 5;
plv_data = [];
BW_freq = [1 3 4 8 8 13 13 30 13 17 20 30 30 50];
for i = 1:30/epoch
    plv_data(:,:,i) = data(:,1+fs*(i-1)*epoch:fs*i*epoch); 
end
for i = 1:length(BW)  
    filtSpec.order = 10;
    filtSpec.range = [BW_freq(2*i-1) BW_freq(2*i)];
    temp_plv = [];
    for j = 1:size(plv_data,1)/2     
        plv = pn_eegPLV(plv_data(2*j-1:2*j,:,:),fs,filtSpec);
        temp_plv = [temp_plv; squeeze(plv(:,1,2,:))'];
    end
    plv_result(:,:,i) = temp_plv(:,:);
end

for bw = 1:length(BW)    
    figure(figure_num)
    figure_num = figure_num + 1;
    for i = 1:length(nameSet)
        subplot(2,2,i);
        p = plot(squeeze(plv_result(4*(i-1)+1:4*i,:,bw))');
        p(1).Color = 'red';
        p(2).Color = 'blue';
        p(3).Color = 'green';
        p(4).Color = 'yellow';
        legend(emotionSet);
        title(char(strcat(nameSet(i)," ", BW(bw)," PLV value")));
    end
    csvwrite(char(strcat(BW(bw)," PLV value.csv")),plv_result(:,:,bw));
end

%% corelation exploration
for j = 1:length(emotionSet)
    for i = 1:length(BW)
        temp_1 = [squeeze(PSD_final(2*j-1,:,i)); squeeze(PSD_final(2*j+7,:,i));squeeze(PSD_final(2*j+15,:,i));squeeze(PSD_final(2*j+23,:,i))];
        temp_2 = [squeeze(PSD_final(2*j,:,i)); squeeze(PSD_final(2*j+8,:,i));squeeze(PSD_final(2*j+16,:,i));squeeze(PSD_final(2*j+24,:,i))];
        COR_channel_1(:,:,i,j) = corr(temp_1');
        COR_channel_2(:,:,i,j) = corr(temp_2');
    end
    temp_3 = [squeeze(Alpha_beta_ratio(2*j-1,:)); squeeze(Alpha_beta_ratio(2*j+7,:));squeeze(Alpha_beta_ratio(2*j+15,:));squeeze(Alpha_beta_ratio(2*j+23,:))];
    temp_4 = [squeeze(Alpha_beta_ratio(2*j,:)); squeeze(Alpha_beta_ratio(2*j+8,:));squeeze(Alpha_beta_ratio(2*j+16,:));squeeze(Alpha_beta_ratio(2*j+24,:))];
    COR_channel_1_ratio_1(:,:,j) = corr(temp_3');
    COR_channel_2_ratio_1(:,:,j) = corr(temp_4');
    temp_5 = [squeeze(Low_beta_high_beta_ratio(2*j-1,:)); squeeze(Low_beta_high_beta_ratio(2*j+7,:));squeeze(Low_beta_high_beta_ratio(2*j+15,:));squeeze(Low_beta_high_beta_ratio(2*j+23,:))];
    temp_6 = [squeeze(Low_beta_high_beta_ratio(2*j,:)); squeeze(Low_beta_high_beta_ratio(2*j+8,:));squeeze(Low_beta_high_beta_ratio(2*j+16,:));squeeze(Low_beta_high_beta_ratio(2*j+24,:))];
    COR_channel_1_ratio_2(:,:,j) = corr(temp_5');
    COR_channel_2_ratio_2(:,:,j) = corr(temp_6');
end











