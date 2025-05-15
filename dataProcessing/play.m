clc; clear; close all;
N = 10000000;
x = rand(4, N);
T = abs(x(1, :) - x(2, :)) + abs(x(3, :) - x(4, :));
histogram(T, 100);
