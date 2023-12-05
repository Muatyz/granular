function [fspan, dft] = dftRead(filename)
%%% 1. Read the data in the dtf.csv file
data = readmatrix(filename);
[Length,~] = size(data);
fspan = data(2:Length,1); dft = data(2:Length,2);
end