function Int = intensity(Amp)
% 求时域信号 Amp 的强度
Amp_squared = Amp.^2;
Int = hilbert(hilbert(Amp_squared));

end