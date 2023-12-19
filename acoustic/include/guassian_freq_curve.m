function curve = guassian_freq_curve(freq_center, sigma)
%%% Guassian distribution: 
%%%
%%% $$
%%% F(f) = (\sigma\sqrt{2\pi})^{-1} e^{-(f - f_{c)^2/(2\sigma^{2})}
%%% $$
%%% $\sigma$: standard deviation; $f_{c}$: freq center;

f = (1:1000)*10^3;
F = 1/(sigma*sqrt(2*pi))*exp(-(f - freq_center).^2/(2*sigma^2)).*exp(1i*2*pi*(f-freq_center));

fs = 2000e3;

%%% plot the freq-domain curve
figure()
plot(f/10^3, F),xlabel("Frequency(kHz)"), ylabel("Magnitude(V)")


%%% IFFT the freq-domain to time domain
curve = ifft(F);

figure()
plot(1:1000, curve)
end