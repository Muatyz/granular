clear;clc;

% Read the original data in *.csv

%%% Read the 46mm distance in air *.csv
[tspan_rn_B_46, Amp_rn_B_46, fspan_rn_B_46, DFT_rn_B_46] = csvRead("大容器-46mm-ref-noise.csv");
[tspan_B_46, Amp_B_46, fspan_B_46, DFT_B_46] = csvRead("大容器-46mm-100hz-10vpp-4us.csv");

[tspan_rn_S_46, Amp_rn_S_46, fspan_rn_S_46, DFT_rn_S_46] = csvRead("小容器-46mm-ref-noise.csv");
[tspan_S_46, Amp_S_46, fspan_S_46, DFT_S_46] = csvRead("小容器-46mm-100hz-10vpp-4us.csv");

%%% Read the 46mm distance in granular media under different stress *.csv
[tspan_S_46_2us_0F, Amp_S_46_2us_0F, fspan_S_46_2us_0F, DFT_S_46_2us_0F] = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-2us-0应力.csv");
[tspan_S_46_2us_F, Amp_S_46_2us_F, fspan_S_46_2us_F, DFT_S_46_2us_F]     = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-2us-应力.csv");

[tspan_S_46_4us_0F, Amp_S_46_4us_0F, fspan_S_46_4us_0F, DFT_S_46_4us_0F] = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-4us-0应力.csv");
[tspan_S_46_4us_F, Amp_S_46_4us_F, fspan_S_46_4us_F, DFT_S_46_4us_F]     = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-4us-应力.csv");

[tspan_S_rn_46_0F, Amp_S_rn_46_0F, fspan_S_rn_46_0F, DFT_S_rn_46_0F] = csvRead("小容器-2mm钢珠-46mm-ref-noise-0应力.csv");
[tspan_S_rn_46_F, Amp_S_rn_46_F, fspan_S_rn_46_F, DFT_S_rn_46_F] =     csvRead("小容器-2mm钢珠-46mm-ref-noise-应力.csv");

% Plot the time-domain signal in subplot form
%%% air boundary effect
figure(1)
subplot(2,1,1)
plot(tspan_B_46*10^3, Amp_B_46,tspan_rn_B_46*10^3, Amp_rn_B_46),xlabel("Time(ms)"),ylabel("Amplitude(V)")
title("L=46mm 不同容器边界的空气信号传播对比")
xlim([0.6,3])
legend("Received Signal","Reference Noise")
subplot(2,1,2)
plot(tspan_S_46*10^3, Amp_S_46,tspan_rn_S_46*10^3, Amp_rn_S_46),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([1.65,3])
legend("Received Signal","Reference Noise")

%%% air boundary effect
figure(2)
subplot(2,1,1)
plot(fspan_B_46, DFT_B_46,fspan_rn_B_46, DFT_rn_B_46),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
title("L=46mm 不同容器边界的空气信号传播对比(FFT)")
legend("Received Signal","Reference Noise")
subplot(2,1,2)
plot(fspan_S_46, DFT_S_46,fspan_rn_S_46, DFT_rn_S_46),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal","Reference Noise")

%%% stress and self weight effect
%%% plot the time domain signal
figure(3)
subplot(3,1,1)
plot(tspan_S_46*10^3, Amp_S_46, tspan_rn_S_46*10^3, Amp_rn_S_46),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([1.65,2.8])
legend("Received Signal(Air)","Reference Noise")
title("检查自重和应力效应的对比")
subplot(3,1,2)
plot(tspan_S_46_4us_0F*10^3, Amp_S_46_4us_0F,tspan_S_rn_46_0F*10^3, Amp_S_rn_46_0F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([0.06,1.05])
legend("Received Signal(Granular, 0 Stress)","Reference Noise")
subplot(3,1,3)
plot(tspan_S_46_4us_F*10^3, Amp_S_46_4us_F, tspan_S_rn_46_F*10^3, Amp_S_rn_46_F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([2.38,3.2])
legend("Received Signal(Granular, Stress)","Reference Noise")

%%% stress and self weight effect
%%% plot the dft
figure(4)
subplot(3,1,1)
plot(fspan_S_46, DFT_S_46,fspan_rn_S_46, DFT_rn_S_46),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(Air)","Reference Noise")
title("检查自重和应力效应的对比(FFT)")
subplot(3,1,2)
plot(fspan_S_46_4us_0F, DFT_S_46_4us_0F,fspan_S_rn_46_0F, DFT_S_rn_46_0F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(Granular, 0 Stress)","Reference Noise")
subplot(3,1,3)
plot(fspan_S_46_4us_F, DFT_S_46_4us_F,fspan_S_rn_46_F, DFT_S_rn_46_F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(Granular, Stress)","Reference Noise")

%% csvReader function
function [tspan, Amp, fspan_half, DFT_half] = csvRead(filename)
%%% 1. Read the data in the *.csv file
data = readmatrix(filename);
[Length, ~] = size(data);
Amp = data(4:Length, 1);
num = Length - 3;
Freq = num / 5 * 10^3;
% delta_t = 1/Freq;

tspan = (0:num - 1) / Freq;
fspan = (0:num - 1) * (Freq / num);

%%% 2. Draw the time-domain signal
%%% Please finish this step on your own *.m file!

%%% 3. FFT. 
%%% Remember that the sampling rate is 2 MHz, so the time step
%%% equals to 1 us.
%%% 
%%% Also remember to use dftRead() function to compare the FFT 
%%% algorithm.

DFT = abs(fft(Amp));
half_num = floor(num / 2);
fspan_half = fspan(1:half_num)/10^3; % Here the unit is 1kHz
DFT_half = DFT(1:half_num);

end
%% DFT table Reader and Visualization
function [fspan, dft] = dftRead(filename)
%%% 1. Read the data in the dtf.csv file
data = readmatrix(filename);
[Length,~] = size(data);
fspan = data(2:Length,1); dft = data(2:Length,2);
end