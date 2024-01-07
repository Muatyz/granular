function [tspan, Amp_raw] = swaeTxtRead(filename)
data = readmatrix(filename);
fs = 2000e3;
[Length,~] = size(data);

Amp_raw = data;
tspan = (0:(Length - 1))/fs;
end