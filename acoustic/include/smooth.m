function signal_smoothed = smooth(signal, fs, window_size, method)
%%% `signal`: shear-stress for example. 
%%% `fs`: sampling frequency set when recording the signal
%%% `window_size`: number of signal points involved in averaging. Don't be
%%% too large to get an appropriate curve
%%% `method`:   - "line"
%%%             - "quadratic"
N = window_size; l = length(signal);

signal_smoothed = zeros(1,l);
radius = floor(N/2);    %%% (o o) o (o o), window = 5, radius = 2
                        %%% (o o o) o (o o o ), window = 6, radius = 3
index_padding = radius;
index_start = index_padding + 1;
index_end = l - radius; %%% 参与平均的最后一个中心点

if method == "line"
    for i = 1:l
        if i < index_start
            signal_smoothed(i) = 0;
        elseif i>= index_start && i <= index_end
            signal_slice = signal((i-radius):(i+radius));
            signal_smoothed(i) = sum(signal_slice)/length(signal_slice);
        else
            signal_smoothed(i) = signal(i);
        end
    end
elseif method == "quadratic"
    for i = 1:l
        if i < index_start
            signal_smoothed(i) = 0;
        elseif i>= index_start && i <= index_end
            slice = signal((i-radius):(i+radius));
            n = length(slice);
            j = (1:n)-1;
            w = 1 - ((j - i).^2*2/l).^2;
            signal_smoothed(i) = sum(w.*slice)/sum(w);
        else
            signal_smoothed(i) = signal(i);
        end
    end
else
    error("Not valid value for `method`! Please check `help` to use the right parameter")
end

end