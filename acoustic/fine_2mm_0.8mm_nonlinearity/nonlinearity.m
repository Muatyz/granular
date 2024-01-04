clear;clc;
addpath('../include/')
%% Preparation

%%% All freq unit is Hz
num = 40; %% the experiment total number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index_filtered = 4; %%% 由于信噪比过低, 在分析波前时需要被提前滤去
%%%%%%%%%%%%%%%%%%  %%% e.g. index_filtered = 3, 说明在进行相似度计算时, 前 3 个数据被滤去. 第
%%%%%%%%%%%%%%%%%%  %%% (3 + 1) 个数据开始分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 2000e3; T = 5 * 10^(-3); sample_length = fs * T; %% Ideal length. Actually, the sample length could be 9997 or something like this
freq_c = 150e3; freq_halfwidth = 5e3;

%%% Create 0s matrix to restore the vectors

Amp_corrected_m = zeros(num,sample_length);
Amp_harmonic = zeros(3,num - index_filtered); %% record each harmonic wave amplitude

efficiency = 0.633653; %%% the response factor at f_c = 150 kHz
                       %%% 可以运行 function_test.m 查看响应效率函数

%% Read data

%%% Read the original signal and correct them with response function
%%% Observe the time of the wave front in Amp_1
%%% 完成这一步之后, 观察波形的信噪比, 再确定 index_filtered. 若信噪都好则可取 0

for i = 1:num
    %%% Read the .csv
    [tspan_raw, Amp_raw, ~, ~] = csvRead(i + ".csv"); 
    
    %%% Response function correction
    [tspan_corrected, Amp_corrected, ~, ~] = response_correct(Amp_raw, "continuous", 0);

    %%% Restore the corrected amplitude vector
    Amp_corrected_m(i, 1:length(Amp_corrected)) = Amp_corrected;

    %%% Plot the ith signal curve on subplot
    figure(1)
    subplot(num,1,i);
    plot(tspan_corrected*10^3, Amp_corrected)
    % xlabel("T(ms)"), ylabel("Amp.(V)");
    legend(i + "th Corr. Signal");
end

%% Similarity Parameter

%%% 注意: 本次实验中, 作为参考的开头信噪比反而是最糟糕的
%%% 因此需要在绘制所有波列信号后, 选取合适的分析序号起始点
%%% 为什么会有交叉积为负的情况? 之前的交叉积计算取最大值忘记取实部了

%%% calculate the cross-relevant function. Choose the 1st is set to be 1.
time_window = 0.8; % Typical time length of the wave train. The unit is ms
index_start = index_filtered + 1;
index_len = num - index_filtered; % 真正可供分析的数据矢量长度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_wavefront_amp_start = 1.186; %%% 序号为 index_start 的波列估计的波前时间. The unit is ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = zeros(1, index_len); 
coefficients_m = zeros(1, index_len); % 创建记录波前和相关系数的空白矩阵

index(1) = time_wavefront_amp_start/10^3 * fs; 
index_1 = index(1);

coefficients_m(1) = 1; %% of course any signal is completely the same with itself

%%% calculate the **neighbor** cross relevant function to find the
%%% wave-front index

for j = 2: index_len
    [coefficient, index_2] = correlate(Amp_corrected_m(index_filtered + j - 1,:), Amp_corrected_m(index_filtered + j,:), index_1, time_window);
    index(j) = index_2; index_1 = index_2;
    coefficients_m(j) = coefficient;
end

%%% Remember find the wavefront first, or you will get compeletly wrong max
%%% in the signal band-pass filtered

%%% Plot the wave front index in figure(1) to check if the auto finding
%%% works well.
%%% Remember to set the unit to 1ms

for j = 1: index_len
figure(1)
subplot(num ,1, index_filtered + j);
xline(index(j)/fs*10^3,               'r','HandleVisibility','off') %% wave-front
xline(index(j)/fs*10^3 + time_window, 'r','HandleVisibility','off') %% wave-tail
end

%%% Plot the sequence order similarity curve
V_input = (1:num)/4;

figure(2)
plot(V_input(index_start:num), coefficients_m,'-o',"LineWidth",1)
xlabel("V_{input}(Vpp)"), ylabel("Similarity Coefficient")

%% Band-pass filter to get each harmonic wave amplitude

time_window_length = fs * time_window / 10^3;  % 比较窗口的长度
waves_m = zeros(index_len, time_window_length);
wave_tspan = (0: (time_window_length - 1))/fs; % Unit is Sec
for i = 1:index_len
    waves_m(i,:) = Amp_corrected_m(i, index(i):(index(i) + time_window_length - 1));
end

for i = 1: index_len
    %%% band pass around 1omega, 2omega and 3omega
    [tspan_filtered_1, Amp_filtered_1] = bandpass(wave_tspan, waves_m(i,:), freq_c*1, freq_halfwidth);
    [tspan_filtered_2, Amp_filtered_2] = bandpass(wave_tspan, waves_m(i,:), freq_c*2, freq_halfwidth);
    [tspan_filtered_3, Amp_filtered_3] = bandpass(wave_tspan, waves_m(i,:), freq_c*3, freq_halfwidth);

    %%% Find the max value in amp_nomega
    [amp_1omega, index_amp_1omega] = max(real(Amp_filtered_1));
    Amp_harmonic(1, i) = amp_1omega;
    
    [amp_2omega, index_amp_2omega] = max(real(Amp_filtered_2));
    Amp_harmonic(2, i) = amp_2omega;
    
    [amp_3omega, index_amp_3omega] = max(real(Amp_filtered_3));
    Amp_harmonic(3, i) = amp_3omega;
end

%% Plot the V_input versus V_iomega^i and the corresponding linear-fit curves

V_eff = V_input * efficiency; 
% V_eff = V_input;
degree = 1; 
V_window = V_eff(index_start:num); %% **可供分析**的信号对应的输入电压矢量
fit_effi = zeros(3, 2); %% 用于存储拟合直线的斜率和截距的空白矩阵

index_cutoff = 13; %%% 发生不可逆转变的实验编号(不含已滤去信号情况下计数, 即在 V_window - Amp_harmonic 图中观察后在填写)
V_cut_off = (index_filtered + index_cutoff)/4*efficiency;

figure(3)
%%% 1w
subplot(3,1,1)
fit_effi(1,:) = polyfit(V_window(1:index_cutoff), Amp_harmonic(1,1:index_cutoff), degree);
slope = fit_effi(1,1);
plot(V_window, Amp_harmonic(1,:),'-o') 
hold on
plot(V_window, polyval(fit_effi(1,:), V_window),'r--',"LineWidth",1)
xline(V_window(index_cutoff), "b")
hold off
xlabel("V_{input}(V)"), ylabel("V_{1\omega}(V)"),
legend("Fundamental","Fit(k = " + slope + ")")

%%% 2w
subplot(3,1,2)
V_window_square = V_window.^2;
fit_effi(2,:) = polyfit(V_window_square(1:index_cutoff), Amp_harmonic(2,1:index_cutoff), degree);
slope = fit_effi(2,1);
plot(V_window_square, Amp_harmonic(2,:),'-o'),
hold on
plot(V_window_square, polyval(fit_effi(2,:), V_window_square),'r--',"LineWidth",1)
xline(V_window_square(index_cutoff), "b")
hold off
xlabel("V_{input}^{2}(V^{2})"), ylabel("V_{2\omega}(V)"), 
legend("2_{nd}", "Fit(k = " + slope + ")")

%%% 3w
subplot(3,1,3)
plot(V_window.^3, Amp_harmonic(3,:),'-o'),xlabel("V_{input}^{3}(V^{3})"), ylabel("V_{3\omega}(V)"), legend("3_{rd}")

%%% Add the cut-off index in the similarity figure to compare the 
%%% reversibility

figure(2)
xline(V_window(index_cutoff), "r")
