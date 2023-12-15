clear;clc;
addpath('../include/')
%% Preparation

%%% All freq unit is Hz
num = 20; %% the experiment total number
fs = 2000e3; T = 5 * 10^(-3); sample_length = fs * T; %% Ideal length. Actually, the sample length could be 9997 or something like this
freq_c = 90e3; freq_halfwidth = 5e3;

%%% Create 0 matrix to restore the vectors

Amp_corrected_m = zeros(num,sample_length);
Amp_harmonic = zeros(3,num);
Index_max_amp_1omega = zeros(1,num);

efficiency = 0.87178; %% the response factor at f_c = 90 kHz

%% Read data

%%% Read the original signal and correct them with response function
%%% Observe the time of the wave front in Amp_1

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
    plot(tspan_corrected*10^3, Amp_corrected),xlabel("Time(ms)"), ylabel("Amplitude(V)"),legend(i + "th Corrected Signal");

end

%% Similarity Parameter

%%% calculate the cross-relevant function. Choose the 1st is set to be 1.
time_window = 1.25; % The unit is ms
time_wavefront_amp1 = 0.4; % The unit is ms
index = zeros(1,num); coefficients_m = zeros(1,num);

index(1) = time_wavefront_amp1/10^3 * fs; 
index_1 = index(1);

coefficients_m(1) = 1; %% of course any signal is completely the same with itself

%%% calculate the **neighbor** cross relevant function
for j = 2 : num
    [coefficient, index_2] = correlate(Amp_corrected_m(j-1,:), Amp_corrected_m(j,:), index_1, time_window);
    index(j) = index_2; index_1 = index_2;
    coefficients_m(j) = coefficient;
end

figure(2)
plot((1:num)/2, coefficients_m,'-o',"LineWidth",1)
xlabel("V_{input}(Vpp)"), ylabel("Similarity Coefficient")

%%% Remember find the wavefront first, or you will get compeletly wrong max
%%% in the signal band-pass filtered

%% Band-pass filter to get each harmonic wave amplitude

time_window_length = fs * time_window / 10^3;
waves_m = zeros(num, time_window_length);
wave_tspan = (0: (time_window_length - 1))/fs; % Unit is Sec
for i = 1:num
    waves_m(i,:) = Amp_corrected_m(i, index(i):(index(i) + time_window_length - 1));
end

for i = 1:num
    %%% band pass around 1omega, 2omega and 3omega
    [tspan_filtered_1, Amp_filtered_1] = bandpass(wave_tspan, waves_m(i,:), freq_c  , freq_halfwidth);
    [tspan_filtered_2, Amp_filtered_2] = bandpass(wave_tspan, waves_m(i,:), freq_c*2, freq_halfwidth);
    [tspan_filtered_3, Amp_filtered_3] = bandpass(wave_tspan, waves_m(i,:), freq_c*3, freq_halfwidth);

    %%% Find the max value in amp_nomega
    [amp_1omega, index_amp_1omega] = max(abs(Amp_filtered_1));
    Amp_harmonic(1, i) = amp_1omega;
    
    [amp_2omega, index_amp_2omega] = max(abs(Amp_filtered_2));
    Amp_harmonic(2, i) = amp_2omega;
    
    [amp_3omega, index_amp_3omega] = max(abs(Amp_filtered_3));
    Amp_harmonic(3, i) = amp_3omega;
end




%%% Plot the V_input versus V_iomega^i and the corresponding linear-fit
%%% curves
figure(3)
V_input = (1:num) /2 * efficiency; degree = 1; 
coefficients_m = zeros(3, 2);
index_cutoff = 13;

subplot(3,1,1)
coefficients_m(1,:) = polyfit(V_input(1:index_cutoff), Amp_harmonic(1,1:index_cutoff), degree);
slope = coefficients_m(1,1);
plot(V_input, Amp_harmonic(1,:),'-o') 
hold on
plot(V_input, polyval(coefficients_m(1,:), V_input),'r--',"LineWidth",1)
xline(V_input(index_cutoff), "b")
hold off
xlabel("V_{input}(V)"), ylabel("V_{1\omega}(V)"),
legend("Fundamental","Fit(k = " + slope + ")")

subplot(3,1,2)
V_input_square = V_input.^2;
coefficients_m(2,:) = polyfit(V_input_square(1:index_cutoff), Amp_harmonic(2,1:index_cutoff), degree);
slope = coefficients_m(2,1);
plot(V_input_square, Amp_harmonic(2,:),'-o'),
hold on
plot(V_input_square, polyval(coefficients_m(2,:), V_input.^2),'r--',"LineWidth",1)
xline(V_input_square(index_cutoff), "b")
hold off
xlabel("V_{input}^{2}(V^{2})"), ylabel("V_{2\omega}(V)"), 
legend("2_{nd}", "Fit(k = " + slope + ")")

subplot(3,1,3)
plot(V_input.^3, Amp_harmonic(3,:),'-o'),xlabel("V_{input}^{3}(V^{3})"), ylabel("V_{3\omega}(V)"), legend("3_{rd}")

%%% Add the cut-off index in the similarity figure to compare the 
%%% reversibility

figure(2)
xline(V_input(index_cutoff), "r")
%% Test the band-pass filter