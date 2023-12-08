function [tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw)
%%% if i = 0, function will plot the response factor function

response_data = readmatrix("../include/response/response.xlsx");
[Length, ~] = size(response_data);

fspan_response = response_data(1:Length,1);
response_amp = response_data(1:Length,2);
signal_amp = response_data(2,7)/2; %% Vpp/2 = Amp

response_sample = sqrt(response_amp./signal_amp);

%%% 扩展响应采样函数到两倍长的频域
fspan_response_double = 1:(2*Length); 
response_sample_double = fspan_response_double;
response_sample_double(1:Length) = response_sample;
for j = 1 : Length
    response_sample_double(Length + j) = response_sample_double(Length - j + 1);
end



%%% Plot the original response factor function on [1,1000] kHz freq domain
figure()
plot(fspan_response, response_sample), yline(1,"r"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
legend("response factor(sample)","Y=1")


%%% FFT
%%% the signal processed in CsvRead() has been cut off so do it
%%% independently
num = length(Amp_raw);
fs = num/(5*10^(-3));

fspan_raw = (0:num - 1) * (fs/num) /10^3; % unit is set to be 1kHz
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