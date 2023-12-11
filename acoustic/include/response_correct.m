function [tspan_corrected, Amp_corrected, fspan_corrected, DFT_corrected] = response_correct(Amp_raw, response_type)

%%% FFT
%%% the signal processed in CsvRead() has been cut off so do it
%%% independently
num = length(Amp_raw);
fs = num/(5*10^(-3));

fspan_raw = (0:num - 1) * (fs/num);
dft_raw = fft(Amp_raw);

if response_type == "continuous"
    response_data = readmatrix("../include/response/continuous_method/response.xlsx",'Range', 'A1:G1009');
    [response_Length, ~] = size(response_data);
    
    fspan_response = (response_data(1:response_Length,1))';
    response_amp = (response_data(1:response_Length,2))';
    signal_amp = 0.5;
    % signal_amp = response_data(2,7)/2; %% Vpp/2 = Amp

    response_sample = sqrt(response_amp./signal_amp);

    %%% 扩展响应采样函数到两倍长的频域
    fspan_response_double = 1:(2*response_Length);
    fspan_response_double(1:response_Length) = fspan_response;
    response_sample_double = fspan_response_double;
    response_sample_double(1:response_Length) = response_sample;
    for j = 1 : response_Length
        fspan_response_double(response_Length + j) = 100 + 2000e3 - fspan_response_double(response_Length - j + 1);
        response_sample_double(response_Length + j) = response_sample_double(response_Length - j + 1);
    end

    %%% Plot the original response factor function on [1,1000] kHz freq domain
    figure()
    plot(fspan_response/10^3, response_sample), yline(1,"r"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
    legend("response factor(sample)","Y=1")


    %%% Interpolation the response factor function according to 'fspan_raw'
    %%% Remember that the fspan_raw is 1-2000khz
    response_inter = interp1(fspan_response_double, response_sample_double, fspan_raw,"spline");
    response_inter(1) = response_inter(2);  %%% fix the nan point

    %%% Correct the freq domain signal with response factor
    dft_corrected = dft_raw./ response_inter;

    %%% Freq domain to Time domain
    Amp_ifft = ifft(dft_corrected);
    Amp_corrected = Amp_ifft;
    tspan_corrected = (0:length(Amp_ifft) - 1)/fs;
    
    %%% Return the required quant
    half_num = floor(num / 2);
    
    fspan_corrected = fspan_raw(1:half_num);
    DFT_corrected = dft_corrected(1:half_num);

elseif response_type == "pulse"
    [tspan_generator, Amp_generator,fspan_generator, dft_generator] = csvRead("../include/response/pulse_method/generator.csv");
    [tspan_received, Amp_received,fspan_received, dft_received] = csvRead("../include/response/pulse_method/received.csv");
    
    %%% Plot the original signal of generator and the received signal
    figure()
    subplot(2,1,1)
    plot(tspan_generator*10^3,Amp_generator),
    xlabel("Time(ms)"),ylabel("Amplitude(V)"),
    legend("Generator Signal"),xlim([0.5,1.5])
    subplot(2,1,2)
    plot(tspan_received*10^3,Amp_received),
    xlabel("Time(ms)"),ylabel("Amplitude(V)"),
    legend("Received Signal"),xlim([0.5,3])

    response_sample = sqrt(abs(dft_received)./abs(dft_generator));

    %%% Plot the response sample function
    figure()
    plot(fspan_generator/10^3,response_sample), yline(1,"r"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
    legend("response factor(sample)","Y=1")


elseif response_type == "noise"
    [tspan_generator, Amp_generator,fspan_generator, dft_generator] = csvRead("../include/response/noise_method/generator.csv");
    [tspan_received, Amp_received,fspan_received, dft_received] = csvRead("../include/response/noise_method/received.csv");

    %%% Plot the original signal of generator and the received signal
    figure()
    subplot(2,1,1)
    plot(tspan_generator*10^3,Amp_generator),
    xlabel("Time(ms)"),ylabel("Amplitude(V)"),
    legend("Generator Signal")
    subplot(2,1,2)
    plot(tspan_received*10^3,Amp_received),
    xlabel("Time(ms)"),ylabel("Amplitude(V)"),
    legend("Received Signal")

    response_sample = sqrt(abs(dft_received)./abs(dft_generator));

    %%% Plot the response sample function
    figure()
    plot(fspan_generator/10^3,response_sample), yline(1,"r"),xlabel("Frequency(kHz)"), ylabel("Magnitude")
    legend("response factor(sample)","Y=1")

else
    error('Not valid parameter in response_correct()!')

end

end