clear;clc;
addpath('../include/')

[tspan, Amp, fspan, dft] = csvRead("小容器-2mm钢珠-36mm-100hz-10vpp-2us-0应力.csv");
Int = intensity(Amp)/50;


figure(1)
plot(tspan*10^3, Amp),xlabel("Time(ms)"),ylabel("Amplitude(V)")

figure(2)
plot(fspan, dft),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")

figure(3)
plot(tspan*10^3, Int),xlabel("Time(ms)"),ylabel("Power(W)")

figure(4)
plot(tspan*10^3, envelope(log10(Int))),xlabel("Time(ms)"),ylabel("log I(a.u.)")
xlim([1.5, 4])

figure(5)
subplot(4,1,1)
plot(tspan*10^3, Amp),xlabel("Time(ms)"),ylabel("Amplitude(V)")
subplot(4,1,2)
plot(fspan, dft),xlabel("Frequency(kHz)"),ylabel("Magnitude(V)")
subplot(4,1,3)
plot(tspan*10^3, Int),xlabel("Time(ms)"),ylabel("Power(W)")
subplot(4,1,4)
plot(tspan*10^3, envelope(log10(Int))),xlabel("Time(ms)"),ylabel("log I(a.u.)")
xlim([1.5, 4])
