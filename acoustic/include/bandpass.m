function [tspan_filtered, Amp_filtered] = bandpass(tspan_in, Amp_in, freq_c, freq_halfwidth)
%%% Advice: every calculation in step should take the standard unit, e.g.
%%% Hz
freq_low = freq_c - freq_halfwidth; freq_high = freq_c + freq_halfwidth;

% 设计带通滤波器
fs = 2000e3;
d = designfilt('bandpassiir', 'FilterOrder', 20, ...
    'HalfPowerFrequency1', freq_low, 'HalfPowerFrequency2', freq_high, ...
    'SampleRate', fs);

% 应用滤波器到信号 S
Amp_filtered = filtfilt(d, Amp_in);
tspan_filtered = tspan_in;
end