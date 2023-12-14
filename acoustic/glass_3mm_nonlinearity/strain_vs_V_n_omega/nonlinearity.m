clear;clc;
addpath('../include/')

%%% All freq unit is Hz
num = 10; %% the experiment total number
fs = 2000e3; T = 5 * 10^(-3); sample_length = fs * T; %% Ideal length. Actually, the sample length could be 9997 or something like this
freq_c = 100e3; freq_halfwidth = 5e3;


%%% Create 0 matrix to restore the vectors

Amp_corrected_m = zeros(num,sample_length);
Amp_harmonic = zeros(3,num);
Index_max_amp_1omega = zeros(1:num);

efficiency = 0.789937; %% the response factor at f_c = 100 kHz

%%% Read the original signal and correct them with response function

%%% Observe the exact index of the wave front


for i = 1:num
    
    [tspan_raw, Amp_raw, ~, ~] = csvRead(i + ".csv"); %% Read the .csv
    
    %%% Response function correction
    [tspan_corrected, Amp_corrected, ~, ~] = response_correct(Amp_raw,"continuous", 0);

    %%% Restore the corrected amplitude vector
    Amp_corrected_m(i, 1:length(Amp_corrected)) = Amp_corrected;

    %%% Plot the ith signal curve on subplot
    figure(1)
    subplot(num,1,i);
    plot(tspan_corrected/10^3, Amp_corrected),xlabel("Time(ms)"), ylabel("Amplitude(V)"),legend(i + "th Corrected Signal");

    %%% band pass around 1omega, 2omega and 3omega
    [tspan_filtered_1, Amp_filtered_1] = bandpass(tspan_corrected, Amp_corrected, freq_c, freq_halfwidth);
    [tspan_filtered_2, Amp_filtered_2] = bandpass(tspan_corrected, Amp_corrected, freq_c*2, freq_halfwidth);
    [tspan_filtered_3, Amp_filtered_3] = bandpass(tspan_corrected, Amp_corrected, freq_c*3, freq_halfwidth);

    %%% Find the max value index in amp_1omega, and look for value of the
    %%% corresponding index in amp_2omega and amp_3omega
    [amp_1omega, index_amp_1omega] = max(abs(Amp_filtered_1));
    Amp_harmonic(1, i) = amp_1omega;
    Index_max_amp_1omega(i) = index_amp_1omega;
    
    Amp_harmonic(2, i) = abs(Amp_filtered_2(index_amp_1omega));
    Amp_harmonic(3, i) = abs(Amp_filtered_3(index_amp_1omega));
end

%%% Plot the V_input versus V_iomega^i
figure(2)
V_input = (1:num) * efficiency;
subplot(3,1,1)
plot(V_input, Amp_harmonic(1,:)),xlabel("V_{input}(V)"), ylabel("V_{1\omega}(V)"), legend("Fundamental")

subplot(3,1,2)
plot(V_input.^2, Amp_harmonic(2,:)),xlabel("V_{input}^{2}(V^{2})"), ylabel("V_{2\omega}(V)"), legend("2_{nd}")

subplot(3,1,3)
plot(V_input.^3, Amp_harmonic(3,:)),xlabel("V_{input}^{3}(V^{3})"), ylabel("V_{3\omega}(V)"), legend("3_{rd}")

%%% calculate the cross-relevant function. Choose the 1st is set to be 1.