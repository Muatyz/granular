clear;clc;
addpath('../include/')

%%% All freq unit is Hz

[tspan_raw, Amp_raw, fspan_raw, dft_raw] = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-2us-0应力.csv");
% [tspan_raw, Amp_raw, fspan_raw, dft_raw] = csvRead("5.csv");

%%% Response function correction
[tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw,"continuous", 1);


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

%%% Freq domain comparisonx
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
response_compare()

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

%%

time_domain_window = 128;
overlap = time_domain_window/2;
nfft = 2 * time_domain_window;

[Spectrogram, f_grid, t_grid]=spectrogram(Amp_corrected, time_domain_window, overlap, nfft, fs, 'yaxis');
half_index = 1:length(f_grid)/2;
f_grid = f_grid(half_index)';
Spectrogram = Spectrogram(half_index, :);

figure()
subplot(2,1,1)
title("Corrected Signal")
plot(tspan_corrected*10^3, Amp_corrected),xlabel("Time(ms)"),ylabel("Amp.(V)")

subplot(2,1,2)
imagesc(t_grid*10^3,f_grid/10^3,abs(Spectrogram));
colormap("jet");colorbar;
title("Freq-Time Heat Map")
xlabel("Time(ms)")
ylabel("Freq(kHz)")

%% 
freq_time_heat_map(tspan_corrected, Amp_corrected,fs)


%% Whole Time File Read
[tspan_all, Amp_all_raw] = csvReadAll("RAEM1_240104211430642.csv");

[tspan_all_c, Amp_all_c, fspan_c, DFT_c] = response_correct(Amp_all_raw,"continuous", 1);

figure()
subplot(2,1,1)
plot(tspan_all_c*10^3, Amp_all_c)
xlabel("Time(ms)"),ylabel("Amplitude");
subplot(2,1,2)
plot(fspan_c/10^3, abs(DFT_c))
xlabel("Freq(kHz)"),ylabel("Magnitude(a.u.)")

%% 调制高斯型的正弦波包

f = 150e3; %% the unit is Hz
n_cycle = 10; 
fs = 2000e3;
delta_t = 1/fs;

%%% Define the parameters of the guassian function
mu = n_cycle/2/f;
sigma = 0.00001; 
%%% generate the curve
t = 0:delta_t:(n_cycle/f);
l = length(t);
gs = 1/(sigma*sqrt(2*pi)) * exp(-(t-mu).^2/(2*sigma^2));
s = sin(2*pi*f*t) .* gs;

fspan = (0:(l - 1)) * (fs/l);
dft = (fft(s))*2;
half_l = floor(l/2);
fspan_half = fspan(1:half_l);
dft_half = dft(1:half_l);

%%% plot the curve
figure()
subplot(2,1,1)
plot(t*10^3, s,t*10^3, gs, t*10^3, envelope(abs(s),10,'peak'),'-',"LineWidth",1)
xlabel("Time(ms)"),ylabel("Magnitude(a.u.)")
legend("\mu="+mu+",\sigma="+sigma,"Guassian Function","Envelope")

subplot(2,1,2)
plot(fspan_half/10^3, abs(dft_half))
xlabel("Frequency(kHz)"),ylabel("Magnitude(a.u.)")
legend("Modulated Siganl DFT")

%% Stress-Strain curve
sampling_rate_force = 0.0005; %% The unit is second.
shearing_rate = 10;
fs_force = 1/sampling_rate_force;
stress = forceRead("stress-strain.csv", fs_force, 0);
num = length(stress);

tspan_f = (0:(num - 1)) / fs_force;

%%% weight-average method to smooth the data
ss_line = smooth(stress, fs_force, 6, "line");
ss_quadratic = smooth(stress, fs_force, 6, "quadratic");

figure()
subplot(3,1,1)
plot(tspan_f, stress)
xlabel("Time(s)"),ylabel("Shear Stress(N)");
legend("Raw Signal at Shearing Rate = "+ shearing_rate +" pps")

subplot(3,1,2)
plot(tspan_f, ss_line)
xlabel("Time(s)"),ylabel("Shear Stress(N)");
legend("Line-Smoothed Signal at Shearing Rate = "+ shearing_rate +" pps")

subplot(3,1,3)
plot(tspan_f, ss_quadratic)
xlabel("Time(s)"),ylabel("Shear Stress(N)");
legend("Quadratic-Smoothed Signal at Shearing Rate = "+ shearing_rate +" pps")
%%% Plot the two curves on the same figure and Zoom it 

figure()
plot(tspan_f, stress, tspan_f, ss_line,tspan_f,ss_quadratic,'--')
xlabel("Time(s)"),ylabel("Shear Stress(N)");
legend("Raw", "Line Averaged", "Quadratic Averaged")

%%% Add a zoomed zone
%%% note: `Image Processing Toolbox` required
zp = BaseZoom();
zp.plot;