function Int = intensity(Amp)
Amp_squared = Amp.^2;
Int = envelope(Amp_squared);
end