function [tspan, Amp, fspan_half, DFT_half] = csvRead(filename)
%%% 1. Read the data in the *.csv file
data = readmatrix(filename);
[Length, ~] = size(data);
Amp = data(4:Length, 1);
num = length(Amp);
Freq = num / 5 * 10^3; % Sampling time is set to be 5 ms
                       % more acceptable way is to set 2000e3 as sampling
                       % freq
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

DFT = (fft(Amp))*2;
half_num = floor(num / 2);
fspan_half = fspan(1:half_num);
DFT_half = DFT(1:half_num);
end