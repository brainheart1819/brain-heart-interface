%% function used for generating file for machine learning
output_data = [];
for i = 1:18
    temp_out_data = squeeze(PSD_final2(1:3,i,:));%./PSD_Mean(1:8,3:6);
    output_data = [output_data; temp_out_data];
    temp_out_data2 = squeeze(PSD_final2(4:6,i,:));%./PSD_Mean(1:8,3:6);
    output_data = [output_data; temp_out_data2];
end

fid = fopen('train_3_26.txt','wt');
[m,n] = size(output_data);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%g,',output_data(i,j));
    end
    
    if mod(i,3) == 1
        fprintf(fid,'%s\n','happy');
    elseif mod(i,3) == 2
        fprintf(fid,'%s\n','sad');
    else
        fprintf(fid,'%s\n','tender');
    end
end