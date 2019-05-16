%% EEG denoising
function re_data = reconstruct_data(filt_data)
[C,L] = wavedec(filt_data,6,'haar');

cA6=appcoef(C,L,'haar',6);
cD6=detcoef(C,L,6);
cD5=detcoef(C,L,5);
cD4=detcoef(C,L,4);
cD3=detcoef(C,L,3);
cD2=detcoef(C,L,2);
cD1=detcoef(C,L,1);

constant = 0.6745;
noise_var_a6 = median(abs(cA6))/constant;
noise_var_d6 = median(abs(cD6))/constant;
noise_var_d5 = median(abs(cD5))/constant;
noise_var_d4 = median(abs(cD4))/constant;
noise_var_d3 = median(abs(cD3))/constant;
noise_var_d2 = median(abs(cD2))/constant;
noise_var_d1 = median(abs(cD1))/constant;
len_a6 = length(cA6);
len_d6 = length(cD6);
len_d5 = length(cD5);
len_d4 = length(cD4);
len_d3 = length(cD3);
len_d2 = length(cD2);
len_d1 = length(cD1);

K1 = 1;
K2 = 3;
K3 = 1.5;

T_a6 = K1 * noise_var_a6 * sqrt(2*log(len_a6));
T_d6 = K1 * noise_var_d6 * sqrt(2*log(len_d6));
T_d5 = K3 * noise_var_d5 * sqrt(2*log(len_d5));
T_d4 = K3 * noise_var_d4 * sqrt(2*log(len_d4));
T_d3 = K1 * noise_var_d3 * sqrt(2*log(len_d3));
T_d2 = K1 * noise_var_d2 * sqrt(2*log(len_d2));
T_d1 = K1 * noise_var_d1 * sqrt(2*log(len_d1));

for i = 1:len_a6
    if abs(cA6(i)) > T_a6
        cA6(i) = T_a6*T_a6/cA6(i);
    end
end

for i = 1:len_d6
    if abs(cD6(i)) > T_d6
        cD6(i) = T_d6*T_d6/cD6(i);
    end
end

for i = 1:len_d5
    if abs(cD5(i)) > T_d5
        cD5(i) = T_d5*T_d5/cD5(i);
    end
end

for i = 1:len_d4
    if abs(cD4(i)) > T_d4
        cD4(i) = T_d4*T_d4/cD4(i);
    end
end

for i = 1:len_d3
    if abs(cD3(i)) > T_d3
        cD3(i) = T_d3*T_d3/cD3(i);
    end
end

for i = 1:len_d2
    if abs(cD2(i)) > T_d2
        cD2(i) = T_d2*T_d2/cD2(i);
    end
end

for i = 1:len_d1
    if abs(cD1(i)) > T_d1
        cD1(i) = T_d1*T_d1/cD1(i);
    end
end

C = [cA6,cD6,cD5,cD4,cD3,cD2,cD1];
re_data = waverec(C,L,'haar');