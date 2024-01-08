function [tspan, Amp_raw] = swaeTxtRead()
filename_list = dir('CH01_Data_*.txt');
data = [];
list_len = length(filename_list);

for i = 1:list_len
    file_name = "CH01_Data_"+num2str(i)+".txt";
    data_slice = (readmatrix(file_name))';
    data = [data, data_slice];
end

fs = 2000e3;
[Length,~] = size(data);

Amp_raw = data;
tspan = (0:(Length - 1))/fs;
end