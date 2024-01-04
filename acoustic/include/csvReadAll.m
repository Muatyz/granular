function [tspan, Amp] = csvReadAll(filename)
%%% 1. Read the data in RAEM1_24XXXXX.csv
data = readmatrix(filename);
[Length, ~] = size(data);
Amp = (data(2:Length, 2))';
num = length(Amp);
fs = 2000e3;
tspan = (0:(num - 1)) / fs;

%%% Plot the Raw Sample curve
figure()
plot(tspan*10^3, Amp)
title("Acoustic Signal in the Whole Run")
xlabel("Time(ms)")
ylabel("Amplitude(V)")
end