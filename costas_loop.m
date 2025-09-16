% clear
% clc

fs = 40e6;
f_in = 1e6+10;
f_d = 4e3;

f_lo_b = 1e6;

N_d = 1000;
N = N_d * round(fs / f_d);

d = 2*(randi([0, 1], [1, N_d]) - 0.5);
% d = ones(1, N_d);
d = repmat(d, floor(fs/f_d), 1);
d = reshape(d, 1, N);

theta0 = pi/6;

kesi = 1/sqrt(2);
bl = 20;

k1 = 2*kesi*(2*bl/(kesi + 1/(4*kesi)));
k2 = (2*bl/(kesi + 1/(4*kesi)))^2;

Tu = 1e-6;
% Tu = 1 / fs;
Dump_N = floor(fs * Tu);

steps = floor(N/Dump_N);
% phase_d = zeros(1, steps);
% phase_d_lpf = zeros(1, steps);
% f_lo_lf = zeros(1, steps);
% f_lo = zeros(1, steps+1);
f_lo(1) = f_lo_b;
p = 0;
s = 0;
% x_i = zeros(1, steps);
% x_q = zeros(1, steps);
% y_i = zeros(1, steps);
% y_q = zeros(1, steps);
% x = zeros(1, steps);
% 
% phi_i = zeros(1, steps);
% phi_lo = zeros(1, steps);

%%
for i = 1:steps
    phi_i(i) = 2*pi*f_in*(i-1)/fs + theta0;
    
    if i == 1
        phi_lo(i) = 0;
    else
        phi_lo(i) = phi_lo(i-1) + 2*pi*f_lo(i)/fs;
    end
    
    phase_d(i) = phi_i(i) - phi_lo(i);
    p = phase_d(i) * k1;
    s = phase_d(i) * Tu * k2 + s;
    f_lo_lf(i) = p + s;
    f_lo(i+1) = f_lo_b + f_lo_lf(i);
end

%%

for i = 1:steps
    if i == 1
        phi_i((i-1)*Dump_N+1:i*Dump_N) = 2*pi*f_in*[(i-1)*Dump_N:i*Dump_N-1]/fs + theta0;
        phi_lo((i-1)*Dump_N+1:i*Dump_N) = 2*pi*f_lo(i)*[(i-1)*Dump_N:i*Dump_N-1]/fs;
    else
        phi_i((i-1)*Dump_N+1:i*Dump_N) = phi_i((i-1)*Dump_N) + 2*pi*f_in*[1:Dump_N]/fs;
        phi_lo((i-1)*Dump_N+1:i*Dump_N) = phi_lo((i-1)*Dump_N) + 2*pi*f_lo(i)*[1:Dump_N]/fs;
    end

    x((i-1)*Dump_N+1:i*Dump_N) = d((i-1)*Dump_N+1:i*Dump_N) .* sin(phi_i((i-1)*Dump_N+1:i*Dump_N));

    x_i((i-1)*Dump_N+1:i*Dump_N) = x((i-1)*Dump_N+1:i*Dump_N) .* sin(phi_lo((i-1)*Dump_N+1:i*Dump_N));
    x_q((i-1)*Dump_N+1:i*Dump_N) = x((i-1)*Dump_N+1:i*Dump_N) .* cos(phi_lo((i-1)*Dump_N+1:i*Dump_N));

    % if i < length(Num)
    %     y_i(i) = Num*[zeros(1, length(Num)-i), x_i(1:i)]' / sum(Num);
    %     y_q(i) = Num*[zeros(1, length(Num)-i), x_q(1:i)]' / sum(Num);
    % else
    %     y_i(i) = Num*x_i(i-length(Num) + 1:i)' / sum(Num);
    %     y_q(i) = Num*x_q(i-length(Num) + 1:i)' / sum(Num);
    % end
    % 
    % x_i_dump(i) = sum(y_i(i)) / Dump_N;
    % x_q_dump(i) = sum(y_q(i)) / Dump_N;

    x_i_dump(i) = sum(x_i((i-1)*Dump_N+1:i*Dump_N)) / Dump_N;
    x_q_dump(i) = sum(x_q((i-1)*Dump_N+1:i*Dump_N)) / Dump_N;

    % phase_d(i) = 1/sqrt(2) * (sign(x_i_dump)*x_q_dump - sign(x_q_dump)*x_i_dump);
    % phase_d(i) = -sign(x_i_dump)*x_q_dump;
    phase_d(i) = x_i_dump(i) * x_q_dump(i) * 4;
    % phase_d(i) = x_i_dump * x_q_dump;
    % if i < length(Num1)
    %     phase_d_lpf(i) = Num1*[zeros(1, length(Num1)-i), phase_d(1:i)]' / sum(Num1);
    % else
    %     phase_d_lpf(i) = Num1*phase_d(i-length(Num1) + 1:i)' / sum(Num1);
    % end
    
    % p = phase_d_lpf(i) * k1;
    % s = phase_d_lpf(i) * Tu * k2 + s;
    p = phase_d(i) * k1;
    s = phase_d(i) * Tu * k2 + s;
    f_lo_lf(i) = p + s;
    f_lo(i+1) = f_lo_b + f_lo_lf(i);
    % f_lo(i+1) = f_lo_b;

end

%%
figure; plot(phase_d);
figure; plot(phi_i - phi_lo)
figure; plot(f_lo_lf);
figure; plot(f_lo);