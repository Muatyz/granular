function [coefficient, index_2] = correlate(Amp_1, Amp_2, index_1, time_window)
%%% Definition of correlation_coefficient
%%% let's say there's 2 signals $S_{i}$ and $S_{j}$, and we want to 
%%% calculate the correlation between them. 
%%% So $ R(i,j) = \frac{C_{i,j}}{\sqrt{C_{i,i}\cdot C_{j,j}}}, 
%%% C_{i,j}(t) = \int_{0}^{T}S_{i}(t)\cdot S_{j}(t)\mathrm{d}t $

%%% Input 2 time domain signal **Amp_1** & **Amp_2**. 
%%%
%%% **index_1** is defined as index of time of the wave_front **Amp_1**, which should 
%%% be observed at the first place. Please calculate it at the first place.
%%%
%%% **time_window** is defined as the typical wave train length, which should 
%%% be observed at the first place. A typical time_window is about 1ms.
%%% (*The unit is ms*, and this function use 2000e3 as sample freq to
%%% calculate something about index
%%% 
%%% Remember to reuse the value index_2 to compare the potential Amp_3 with
%%% Amp_2
%%%
%%% Output the tau-gamma(tau) every sinle time. Use % to comment the plot()
%%% codes to avoid output too many figures


%%% Traversal the whole time(usually [0, 5ms-time_window]) to find the
%%% index of which the coefficient get maximum

fs = 2000e3;
time_window_length = fs * time_window / 10^3;
wave_1 = Amp_1(index_1 : (index_1 + time_window_length - 1));
index_2_m = 1 : (length(Amp_1) - time_window_length + 1);
coeffi_m = 1 : (length(Amp_1) - time_window_length + 1);

for i = 1 : length(index_2_m)
    wave_2 = Amp_2(i : (i + time_window_length - 1));
    coeffi_m(i) = sum(wave_1 .* wave_2)/ (sqrt(sum(wave_1.^2)) * sqrt(sum(wave_2.^2)));
end

%%% comment this if you dont want to output it every single time
%figure()
%plot(index_2_m/fs *10^3, coeffi_m), xlabel("\tau(ms)"),ylabel("\Gamma(\tau)")

[coefficient, index_2] = max(coeffi_m);