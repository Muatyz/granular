clear;clc;
addpath('../include/')

%%% All freq unit is Hz

[tspan_raw, Amp_raw, fspan_raw, dft_raw] = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-2us-0应力.csv");
% Int = intensity(Amp_raw)/50;

%%% Response function correction
[tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct_dev(Amp_raw);


%%% High-Pass Filter
fs = 2000e3; % Sample freq
fc_highpass = 300e3; % High-Pass cut-off freq

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

figure(6)
%%% Original Signal
subplot(6,1,1)
plot(tspan_raw*10^3, Amp_raw),xlabel("Time(ms)"),ylabel("Amplitude(V)")
subplot(6,1,2)
plot(fspan_raw/10^3,abs(dft_raw)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")

%%% Corrected Signal
subplot(6,1,3)
plot(tspan_corrected*10^3, Amp_corrected),xlabel("Time(ms)"),ylabel("Amplitude(V)")
subplot(6,1,4)
plot(fspan_corrected/10^3, abs(DFT_corrected)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")

%%% Filter the corrected signal
subplot(6,1,5)
plot(tspan_corrected*10^3, Amp_corrected_filtered),xlabel("Time(ms)"),ylabel("Amplitude(V)")
subplot(6,1,6)
plot(fspan_corrected_filtered/10^3, abs(dft_corrected_filtered)),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")

%%
function [tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct_dev(Amp_raw)

response_data = readmatrix("../include/response/response.xlsx");
[Length, ~] = size(response_data);

fspan_response = response_data(1:Length,1);
response_amp = response_data(1:Length,2);
signal_amp = response_data(2,7)/2; %% Vpp/2 = Amp

response_sample = sqrt(response_amp./signal_amp);

%%% 扩展响应采样函数到两倍长的频域
fspan_response_double = (1:(2*Length))*1000; 
response_sample_double = fspan_response_double;
response_sample_double(1:Length) = response_sample;
for j = 1 : Length
    response_sample_double(Length + j) = response_sample_double(Length - j + 1);
end

%%% Plot the original response factor function on [1,1000] kHz freq domain
figure()
plot(fspan_response/10^3, response_sample), yline(1,"r"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
legend("response factor(sample)","Y=1")


%%% FFT
%%% the signal processed in CsvRead() has been cut off so do it
%%% independently
num = length(Amp_raw);
fs = num/(5*10^(-3));

fspan_raw = (0:num - 1) * (fs/num);
dft_raw = fft(Amp_raw);

%%% Interpolation the response factor function according to 'fspan_raw'
%%% Remember that the fspan_raw is 1-2000khz
response_inter = interp1(fspan_response_double, response_sample_double, fspan_raw,"linear");

%%% Correct the freq domain signal with response factor
dft_corrected = dft_raw./ response_inter;

%%% Freq domain to Time domain
Amp_ifft = ifft(dft_corrected);Amp_corrected = Amp_ifft;
tspan_corrected = (0:length(Amp_ifft) - 1)/fs;

%%% Return the required quant
half_num = floor(num / 2);

fspan_corrected = fspan_raw(1:half_num);
DFT_corrected = dft_corrected(1:half_num);
end