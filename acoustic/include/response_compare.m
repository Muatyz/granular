function response_compare()
%%% use the function to compare response **sample** function measured with
%%% different methods

%% continuous sine waves exciting the transducer
response_data = readmatrix("../include/response/continuous_method/response.xlsx",'Range', 'A1:G1009');
[response_Length, ~] = size(response_data);
    
fspan_response = (response_data(1:response_Length,1))';
response_amp = (response_data(1:response_Length,2))';
signal_amp = 0.5;
% signal_amp = response_data(2,7)/2; %% Vpp/2 = Amp

response_sample = sqrt(response_amp./signal_amp);

figure()
plot(fspan_response/10^3, response_sample,'LineWidth',2.5), yline(1,"g"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
hold on
%% pulse waves
[~, ~,fspan_generator, dft_generator] = csvRead("../include/response/pulse_method/generator.csv");
[~, ~,~, dft_received] = csvRead("../include/response/pulse_method/received.csv");

response_sample = sqrt(abs(dft_received)./abs(dft_generator));
plot(fspan_generator/10^3,response_sample),xlabel("Frequency(kHz)"), ylabel("Magnitude"),legend("response factor(Pulse Waves)")
%% theoretical 'white' noise
[~, ~,fspan_generator, dft_generator] = csvRead("../include/response/noise_method/generator.csv");
[~, ~,~, dft_received] = csvRead("../include/response/noise_method/received.csv");

response_sample = sqrt(abs(dft_received)./abs(dft_generator));
plot(fspan_generator/10^3,response_sample),xlabel("Frequency(kHz)"), ylabel("Magnitude"),legend("response factor(Noise Waves)")
xline(400,"b"),ylim([0,7])
legend("response factor(Continuous Sine Wave)","Y=1","response factor(Pulse Waves)","response factor(Noise Waves)","f = 400kHz")
hold off
end