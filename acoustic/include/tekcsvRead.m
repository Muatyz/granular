function [fs, tspan, Amp, arrival_time]=tekcsvRead(filename)
data = readmatrix(filename);
[Length, ~] = size(data);

tspan = (data(:, 4));
Amp = (data(:, 5));

fs = 1/data(2,2);
arrival_time = 0;
end