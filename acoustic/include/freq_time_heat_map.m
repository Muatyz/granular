function freq_time_heat_map(tspan_corrected, Amp_corrected, fs)
time_domain_window = 128;
overlap = time_domain_window/2;
nfft = 256;

[Spectrogram, f_span, t_span]=spectrogram(Amp_corrected, time_domain_window, overlap, nfft, fs, 'yaxis');
half_index = 1:length(f_span)/2;
f_span = f_span(half_index)';
Spectrogram = Spectrogram(half_index, :);

figure()
subplot(2,1,1)
title("Corrected Signal")
plot(tspan_corrected*10^3, Amp_corrected),xlabel("Time(ms)"),ylabel("Amp.(V)")

subplot(2,1,2)
imagesc(t_span*10^3,f_span/10^3,abs(Spectrogram));
colormap("jet");colorbar;
title("Freq-Time Heat Map")
xlabel("Time(ms)")
ylabel("Freq(kHz)")

end