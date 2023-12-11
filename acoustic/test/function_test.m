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
fspan = (0:num - 1) * (Freq / num);
fspan_corrected_filtered = fspan(1:half_num); % Here the unit is 1kHz
dft_corrected_filtered = DFT_corrected_filtered(1:half_num);

% figure(1)
% plot(tspan*10^3, Amp),xlabel("Time(ms)"),ylabel("Amplitude(V)")
% 
% figure(2)
% plot(fspan, dft),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")
% 
% figure(3)
% plot(tspan*10^3, Int),xlabel("Time(ms)"),ylabel("Power(W)")
% 
% figure(4)
% plot(tspan*10^3, envelope(log10(Int))),xlabel("Time(ms)"),ylabel("log I(a.u.)")
% xlim([1.5, 4])
% 
% figure(5)
% subplot(4,1,1)
% plot(tspan*10^3, Amp),xlabel("Time(ms)"),ylabel("Amplitude(V)")
% subplot(4,1,2)
% plot(fspan, dft),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")
% subplot(4,1,3)
% plot(tspan*10^3, Int),xlabel("Time(ms)"),ylabel("Power(W)")
% subplot(4,1,4)
% plot(tspan*10^3, envelope(log10(Int))),xlabel("Time(ms)"),ylabel("log I(a.u.)")
% xlim([1.5, 4])

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

figure(3)
Int_corrected_filtered = intensity(abs(Amp_corrected_filtered))/50;
plot(tspan_raw*10^3, Int_corrected_filtered),xlabel("Time(ms)"),ylabel("Power(W)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)"),xlim([1.7,3])

figure(4)
plot(tspan_raw*10^3, log(Int_corrected_filtered)),xlabel("Time(ms)"),ylabel("log I(a.u.)"),legend("Filtered(>" + fc_highpass/10^3 + "kHz)"),xlim([1.7,3])

%% Test the pulse_method

% [tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw,"pulse");
response_compare()