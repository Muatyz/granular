clear;clc;

% Read the original data in *.csv
[tspan_rn, Amp_rn, fspan_rn, DFT_rn] = csvRead("reference_noise.csv");
[tspan_1, Amp_1, fspan_1, DFT_1] = csvRead("1_cycle.csv");
[tspan_5, Amp_5, fspan_5, DFT_5] = csvRead("5_cycle.csv");
[tspan_10, Amp_10, fspan_10, DFT_10] = csvRead("10_cycle.csv");

% Read the DFT data calculated by SWAE
[fspan_rn_swae, dft_rn_swae] = dftRead("reference_noise_dft.csv");
[fspan_1_swae, DFT_1_swae] = dftRead("1_cycle_dft.csv");
[fspan_5_swae, DFT_5_swae] = dftRead("5_cycle_dft.csv");
[fspan_10_swae, DFT_10_swae] = dftRead("10_cycle_dft.csv");

% Plot the time-domain signal in subplot form
figure(1)
subplot(3,1,1)
plot(tspan_10*10^3, Amp_10,tspan_rn*10^3, Amp_rn),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([2.5,5])
legend("10 Cycle","Reference Noise")
subplot(3,1,2)
plot(tspan_5*10^3, Amp_5,tspan_rn*10^3, Amp_rn),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([1.5,4])
legend("5 Cycle","Reference Noise")
subplot(3,1,3)
plot(tspan_1*10^3, Amp_1,tspan_rn*10^3, Amp_rn),xlabel("Time(ms)"),ylabel("Amplitude(V)")
xlim([2.5,5])
legend("1 Cycle","Reference Noise")

% Plot the freq-domain signal in subplot form
figure(2)
subplot(4,1,1)
plot(fspan_10, DFT_10,fspan_10_swae,DFT_10_swae),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
legend("10 Cycle", "10 Cycle-SWAE")
subplot(4,1,2)
plot(fspan_5, DFT_5,fspan_5_swae,DFT_5_swae),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
legend("5 Cycle", "5 Cycle-SWAE")
subplot(4,1,3)
plot(fspan_1, DFT_1,fspan_1_swae,DFT_1_swae),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
legend("1 Cycle", "1 Cycle-SWAE")
subplot(4,1,4)
plot(fspan_rn, DFT_rn,fspan_rn_swae,dft_rn_swae),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
legend("Reference Noise", "Reference Noise-SWAE")

% Plot the freq-domain signal in subplot form
figure(3)
subplot(4,1,1)
plot(fspan_10, DFT_10./(max(DFT_10)),fspan_10_swae,DFT_10_swae./(max(DFT_10_swae))),xlabel("Frequency(kHz)"),ylabel("Amp%")
legend("10 Cycle", "10 Cycle-SWAE")
subplot(4,1,2)
plot(fspan_5, DFT_5./(max(DFT_5)),fspan_5_swae,DFT_5_swae./(max(DFT_5_swae))),xlabel("Frequency(kHz)"),ylabel("Amp%")
legend("5 Cycle", "5 Cycle-SWAE")
subplot(4,1,3)
plot(fspan_1, DFT_1./(max(DFT_1)),fspan_1_swae,DFT_1_swae./(max(DFT_1_swae))),xlabel("Frequency(kHz)"),ylabel("Amp%")
legend("1 Cycle", "1 Cycle-SWAE")
subplot(4,1,4)
plot(fspan_rn, DFT_rn./max(DFT_rn),fspan_rn_swae,dft_rn_swae./(max(dft_rn_swae))),xlabel("Frequency(kHz)"),ylabel("Amp%")
legend("Reference Noise", "Reference Noise-SWAE")

% % Plot the corrected Response Signal
% figure(4)
% subplot(3,1,1)
% plot(fspan_10, DFT_10-DFT_rn),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
% legend("10 Cycle")
% subplot(3,1,2)
% plot(fspan_5, DFT_5-DFT_rn),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
% legend("5 Cycle")
% subplot(3,1,3)
% plot(fspan_1, DFT_1-DFT_rn),xlabel("Frequency(kHz)"),ylabel("Amp(V)")
% legend("1 Cycle")

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