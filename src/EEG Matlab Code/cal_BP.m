%% The function used for calculate band power  
function BP = cal_BP(data, fs)
delta = [1,3];
theta = [4,8];
alpha = [8,13];
beta = [13,30];
low_beta = [13,17];
high_beta = [20,30];
gamma = [30,50];

index = floor(length(data)/1000);
start = 1;
end_  = 1000; 

DeltaBP = [];
ThetaBP = [];
AlphaBP = [];
BetaBP  = [];
Low_BetaBP = [];
High_BetaBP = [];
GammaBP = [];

for i = 1:index
    delta_band = bandpower(data(start:end_), fs, delta)/2;
    theta_band = bandpower(data(start:end_), fs, theta)/4;
    alpha_band = bandpower(data(start:end_), fs, alpha)/5;
    beta_band  = bandpower(data(start:end_), fs, beta)/17;
    beta_band_low = bandpower(data(start:end_), fs, low_beta)/4;
    beta_band_high = bandpower(data(start:end_), fs, high_beta)/10;
    gamma_band = bandpower(data(start:end_), fs, gamma)/20;
    DeltaBP = [DeltaBP, delta_band];
    ThetaBP = [ThetaBP, theta_band];
    AlphaBP = [AlphaBP, alpha_band];
    BetaBP  = [BetaBP,  beta_band ];
    Low_BetaBP = [Low_BetaBP, beta_band_low];
    High_BetaBP = [High_BetaBP, beta_band_high];
    GammaBP = [GammaBP, gamma_band];
    start = i*1000;
    end_  = end_ + 1000;
end
BP = [DeltaBP; ThetaBP; AlphaBP; BetaBP; Low_BetaBP; High_BetaBP; GammaBP];
