clear;clc;

% Read the original data in *.csv

%%% Read the 46mm distance in granular media under different stress *.csv
[tspan_S_46_4us_0F, Amp_S_46_4us_0F, fspan_S_46_4us_0F, DFT_S_46_4us_0F] = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-4us-0应力.csv");
[tspan_S_46_4us_F, Amp_S_46_4us_F, fspan_S_46_4us_F, DFT_S_46_4us_F]     = csvRead("小容器-2mm钢珠-46mm-100hz-10vpp-4us-应力.csv");

[tspan_S_rn_46_0F, Amp_S_rn_46_0F, fspan_S_rn_46_0F, DFT_S_rn_46_0F] = csvRead("小容器-2mm钢珠-46mm-ref-noise-0应力.csv");
[tspan_S_rn_46_F, Amp_S_rn_46_F, fspan_S_rn_46_F, DFT_S_rn_46_F] =     csvRead("小容器-2mm钢珠-46mm-ref-noise-应力.csv");

%%% Read the 36mm distance in granular media under different stress *.csv
[tspan_S_36_4us_0F, Amp_S_36_4us_0F, fspan_S_36_4us_0F, DFT_S_36_4us_0F] = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-4us-0应力.csv");
[tspan_S_36_4us_F, Amp_S_36_4us_F, fspan_S_36_4us_F, DFT_S_36_4us_F]     = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-4us-应力.csv");

[tspan_S_rn_36_0F, Amp_S_rn_36_0F, fspan_S_rn_36_0F, DFT_S_rn_36_0F] = csvRead("小容器-2mm钢珠-36mm-ref-noise-0应力.csv");
[tspan_S_rn_36_F, Amp_S_rn_36_F, fspan_S_rn_36_F, DFT_S_rn_36_F] =     csvRead("小容器-2mm钢珠-36mm-ref-noise-应力.csv");

% Plot the time-domain signal in subplot form

%%% stress and self weight effect
%%% plot the time domain signal
figure(1)
subplot(2,2,1)
plot(tspan_S_36_4us_0F*10^3, Amp_S_36_4us_0F, tspan_S_rn_36_0F*10^3, Amp_S_rn_36_0F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([1.85,3.2])
legend("Received Signal(0 Stress,36 mm)","Reference Noise")
subplot(2,2,2)
plot(tspan_S_36_4us_F*10^3, Amp_S_36_4us_F, tspan_S_rn_36_F*10^3, Amp_S_rn_36_F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([0.5,1.5])
legend("Received Signal(Stress,36 mm)","Reference Noise")
subplot(2,2,3)
plot(tspan_S_46_4us_0F*10^3, Amp_S_46_4us_0F, tspan_S_rn_46_0F*10^3, Amp_S_rn_46_0F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([0.05,1.05])
legend("Received Signal(0 Stress,46 mm)","Reference Noise")
subplot(2,2,4)
plot(tspan_S_46_4us_F*10^3, Amp_S_46_4us_F, tspan_S_rn_46_F*10^3, Amp_S_rn_46_F),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([2.38,3.48])
legend("Received Signal(Stress,46 mm)","Reference Noise")

%%% stress and self weight effect
%%% plot the dft
figure(2)
subplot(2,2,1)
plot(fspan_S_36_4us_0F,DFT_S_36_4us_0F,fspan_S_rn_36_0F,DFT_S_rn_36_0F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(0 Stress,36 mm)","Reference Noise")
subplot(2,2,2)
plot(fspan_S_36_4us_F,DFT_S_36_4us_F,fspan_S_rn_36_F,DFT_S_rn_36_F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(Stress,36 mm)","Reference Noise")
subplot(2,2,3)
plot(fspan_S_46_4us_0F,DFT_S_46_4us_0F,fspan_S_rn_46_0F,DFT_S_rn_46_0F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(0 Stress,46 mm)","Reference Noise")
subplot(2,2,4)
plot(fspan_S_46_4us_F,DFT_S_46_4us_F,fspan_S_rn_46_F,DFT_S_rn_46_F),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
xlim([50,250])
legend("Received Signal(Stress,46 mm)","Reference Noise")
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