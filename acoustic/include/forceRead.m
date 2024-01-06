function stress = forceRead(filename, fs_force, show_curve)
%%% Recommended sampling frequency of force sensor is 1/0.0005 = 2000 Hz
%%% A typical shearing process lasts for about 1000s, so the data set
%%% would be large.
%%%
%%% 'show_curve': set 1 for plotting the shear stress-(strain)
data = readmatrix(filename);
[Length, ~] = size(data);
stress = (data(15:Length, 1))';
num = length(stress);

tspan = (0:(num - 1)) / fs_force;
%%% Set a time window to smooth the data.

if show_curve == 1
    figure()
    plot(tspan, stress)
    xlabel("Time(s)"),ylabel("Shear Stress(N)");
elseif show_curve == 0
    
else 
    error("Not valid flag value for 'show_curve'! Please check it again.")
end