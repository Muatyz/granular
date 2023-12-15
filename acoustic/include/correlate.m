function coefficient = correlate(tspan_1,Amp_1, tspan_2, Amp_2)
%%% Definition of correlation_coefficient
%%% let's say there's 2 signals $S_{i}$ and $S_{j}$, and we want to 
%%% calculate the correlation between them. 
%%% So $ R(i,j) = \frac{C_{i,j}}{\sqrt{C_{i,i}\cdot C_{j,j}}}, 
%%% C_{i,j}(t) = \int_{0}^{T}S_{i}(t)\cdot S_{j}(t)\mathrm{d}t $

%%% Plot the signal in time domain, calculate the delayed time $\tau$

figure()
plot(tspan_1,Amp_1, tspan_2, Amp_2)

%%% Find the correct delayed time $\tau$

end