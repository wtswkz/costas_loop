clear
clc

fs = 40e6;
f_in = 10e6;
f_d = 4e3;

f_lo_b = 10e6;

N_d = 8000;
N = N_d * round(fs / f_d);

d = 2*(randi([0, 1], [1, N_d]) - 0.5);
d = repmat(d, floor(fs/f_d), 1);
d = reshape(d, 1, N);

theta0 = pi/12;

kesi = 1/sqrt(2);
bl = 20;

k1 = 2*kesi*(2*bl/(kesi + 1/(4*kesi))); % tau2carr / tau1carr
k2 = (2*bl/(kesi + 1/(4*kesi)))^2; % 1 / tau1carr

Tu = 0.001;
Dump_N = floor(fs * Tu);

steps = floor(N/Dump_N);
phase_d = zeros(1, steps);
f_lo_lf = zeros(1, steps);
f_lo = zeros(1, steps+1);
f_lo(1) = f_lo_b;
p = 0;
s = 0;

for i = 1:steps
    x = d((i-1)*Dump_N+1:i*Dump_N) .* cos(2*pi*f_in*[(i-1)*Dump_N+1:i*Dump_N]/fs + theta0);
    % x = cos(2*pi*f_in*[(i-1)*Dump_N:i*Dump_N-1]/fs + theta0);
    x_i = x .* cos(2*pi*f_lo(i)*[(i-1)*Dump_N:i*Dump_N-1]/fs);
    x_q = x .* sin(2*pi*f_lo(i)*[(i-1)*Dump_N:i*Dump_N-1]/fs);
    % x_i = x(1:2:Dump_N);
    % x_q = x(2:2:Dump_N);
    % x_c = (x_i + j*x_q).*(cos(2*pi*f_lo(i)*[(i-1)*Dump_N/2:i*Dump_N/2-1]/fs) - j*sin(2*pi*f_lo(i)*[(i-1)*Dump_N/2:i*Dump_N/2-1]/fs));
    % x_i = d((i-1)*Dump_N+1:i*Dump_N) .* cos(2*pi*(f_in - f_lo(i))*[(i-1)*Dump_N+1:i*Dump_N]/fs + theta0);
    % x_q = d((i-1)*Dump_N+1:i*Dump_N) .* sin(2*pi*(f_in - f_lo(i))*[(i-1)*Dump_N+1:i*Dump_N]/fs + theta0);
    % x_i_dump = sum(real(x_c)) / Dump_N * 2;
    % x_q_dump = sum(imag(x_c)) / Dump_N * 2;
    x_i_dump = sum(x_i) / Dump_N;
    x_q_dump = sum(x_q) / Dump_N;

    % phase_d(i) = sign(x_q_dump)*x_i_dump - sign(x_i_dump)*x_q_dump;
    phase_d(i) = x_i_dump * x_q_dump;
    p = phase_d(i) * k1;
    s = phase_d(i) * Tu * k2 + s;
    f_lo_lf(i) = p + s;
    f_lo(i+1) = f_lo_b + f_lo_lf(i); 

end

%%
figure; plot(phase_d);
figure; plot(f_lo_lf);
figure; plot(f_lo);