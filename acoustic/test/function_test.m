clear;clc;
addpath('../include/')

%%% All freq unit is Hz

[tspan_raw, Amp_raw, fspan_raw, dft_raw] = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-2us-0应力.csv");

%%% Response function correction
[tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw,"continuous");


%%% High-Pass Filter
fs = 2000e3; % Sample freq
fc_highpass = 140e3; % High-Pass cut-off freq

order = 4;
[b,a] = butter(order,fc_highpass/(fs/2), 'high');
Amp_corrected_filtered = filter(b, a, Amp_corrected);

DFT_corrected_filtered = abs(fft(Amp_corrected_filtered));
num = length(Amp_corrected_filtered);
Freq = num / 5 * 10^3;
half_num = floor(num / 2);
fspan = (0:(num - 1)) * (Freq / num);
fspan_corrected_filtered = fspan(1:half_num); % Here the unit is 1kHz
dft_corrected_filtered = DFT_corrected_filtered(1:half_num)*2;

%%% Time domain comparison
figure(1)
subplot(3,1,1)
plot(tspan_raw*10^3, Amp_raw),xlabel("Time(ms)"),ylabel("Amplitude(V)"),legend("Original Signal"),xlim([1.7,3])
subplot(3,1,2)
plot(tspan_corrected*10^3, Amp_corrected),xlabel("Time(ms)"),ylabel("Amplitude(V)"),legend("Corrected Signal"),xlim([1.7,3])
subplot(3,1,3)
plot(tspan_corrected*10^3, Amp_corrected_filtered),xlabel("Time(ms)"),ylabel("Amplitude(V)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)"),xlim([1.7,3])

%%% Freq domain comparison
figure(2)
subplot(3,1,1)
plot(fspan_raw/10^3,abs(dft_raw)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)"),legend("Original Signal")
subplot(3,1,2)
plot(fspan_corrected/10^3, abs(DFT_corrected)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)"),legend("Corrected Signal")
subplot(3,1,3)
plot(fspan_corrected_filtered/10^3, abs(dft_corrected_filtered)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)")


%%% Intensity Plot

%%% I
figure(3)
Int_corrected_filtered = intensity(abs(Amp_corrected_filtered))/50;
plot(tspan_raw*10^3, Int_corrected_filtered),xlabel("Time(ms)"),ylabel("Power(W)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)"),xlim([1.7,3])


%%% extract the time interval signal and fit the line
start_time = 2.788;end_time = 3.716;%% Unit is 1 ms. Plot the log plot first to find the proper time interval
start_index = find(tspan_raw*10^3 >= start_time, 1);
end_index = find(tspan_raw*10^3 <= end_time, 1, 'last');

t_interval = tspan_raw(start_index:end_index)*10^3;
log_Int = log10(Int_corrected_filtered);
log_Int_interval = log_Int(start_index:end_index);

degree = 1;
coefficients = polyfit(t_interval, log_Int_interval, degree);
slope = coefficients(1);

%%% log I
figure(4)
plot(tspan_raw*10^3, log10(Int_corrected_filtered));
hold on
plot(t_interval, polyval(coefficients, t_interval),'r-',"LineWidth",2)
hold off
xlabel("Time(ms)"),ylabel("log10 I(a.u.)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)", "Fitted Line, \tau="+ abs((-log10(exp(1)))/slope) + "ms");
xlim([1.7,4])


%% Test the pulse_method

% [tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw,"pulse");
% response_compare()

%% Test the band pass filter
freq_c = 100e3; freq_halfwidth = 5e3;
freq_low = freq_c - freq_halfwidth; freq_high = freq_c + freq_halfwidth;
[tspan_filtered_1, Amp_filtered_1] = bandpass(tspan_corrected, Amp_corrected, freq_c, freq_halfwidth);
[tspan_filtered_2, Amp_filtered_2] = bandpass(tspan_corrected, Amp_corrected, freq_c*2, freq_halfwidth);
[tspan_filtered_3, Amp_filtered_3] = bandpass(tspan_corrected, Amp_corrected, freq_c*3, freq_halfwidth);


figure(6)
subplot(2,1,1)
plot(tspan_raw*10^3, Amp_raw),xlabel("Time(ms)"),ylabel("Amplitude(V)"),legend("Corrected Signal"),xlim([1.7,3])
subplot(2,1,2)
plot(tspan_filtered_1*10^3, Amp_filtered_1)
hold on
plot(tspan_filtered_1*10^3, envelope(abs(Amp_filtered_1),20,'peak'),'-',"LineWidth",1),xlabel("Time(ms)"),ylabel("Amplitude(V)");
hold off
legend("Band-Pass Filtered Signal(" + freq_low/10^3 + "~" + freq_high/10^3 + " kHz)", "Envelope"),xlim([1.7,3])

figure(7)
subplot(3,1,1)
plot(tspan_filtered_1*10^3, Amp_filtered_1)
hold on
plot(tspan_filtered_1*10^3, envelope(abs(Amp_filtered_1),20,'peak'),'-',"LineWidth",1),xlabel("Time(ms)"),ylabel("Amplitude(V)");
hold off
legend("Band-Pass Filtered Signal(" + freq_low/10^3 + "~" + freq_high/10^3 + " kHz)", "Envelope"),xlim([1.7,3])
subplot(3,1,2)
plot(tspan_filtered_2*10^3, Amp_filtered_2)
hold on
plot(tspan_filtered_2*10^3, envelope(abs(Amp_filtered_2),20,'peak'),'-',"LineWidth",1),xlabel("Time(ms)"),ylabel("Amplitude(V)");
hold off
legend("Band-Pass Filtered Signal(" + (freq_c*2 - freq_halfwidth)/10^3 + "~" + (freq_c*2 + freq_halfwidth)/10^3 + " kHz)", "Envelope"),xlim([1.7,3])
subplot(3,1,3)
plot(tspan_filtered_3*10^3, Amp_filtered_3)
hold on
plot(tspan_filtered_3*10^3, envelope(abs(Amp_filtered_3),20,'peak'),'-',"LineWidth",1),xlabel("Time(ms)"),ylabel("Amplitude(V)");
hold off
legend("Band-Pass Filtered Signal(" + (freq_c*3 - freq_halfwidth)/10^3 + "~" + (freq_c*3 + freq_halfwidth)/10^3 + " kHz)", "Envelope"),xlim([1.7,3])