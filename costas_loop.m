clear
clc

%%
fs = 40e6;
f_in = 10e6;
f_d = 4e3;

f_lo = 10e6;

N_d = 10;
N = N_d * round(fs / f_d);

d = 2*(randi([0, 1], [1, N_d]) - 0.5);
d = repmat(d, round(fs/f_d), 1);
d = reshape(d, 1, N_d * round(fs/f_d));

theta0 = pi/6;
x = d * sin(2*pi*f_in*[0:N-1]/fs + theta0);

