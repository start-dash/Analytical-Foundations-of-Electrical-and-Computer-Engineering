%% Project 2 - Problem 1


%% part c (plotting the magnitude of the transfer function)
omega = logspace(-2,4,10000); %% c.a.

L = 0.01;
R = 3.38;
K = 0.029;
J = 2e-4;
beta = 0.5e-5;

s = j.*omega;

%% vector for H(w)

tau = L/R; %% for an RL circuit
H = K ./ ((R*J*s)+(beta*R)+(L*J*(s).^2)+(L*beta*s)+(K).^2);

figure(1);
semilogx(omega,20*log10(abs(H))); %% c.c.

title('Transfer Function Magnitude');
xlabel('Frequency');
ylabel('Magnitude');


%% part i
sw = 106; % switching frequency
Tsw = 2*pi/sw; % period

delta_t = Tsw/10000;
t = [0:delta_t:10];

Va = 6*square(sw*t, 61.63) + 6;

omega = [];
I_a = [];

for n = 1:1:length(t)-1
    if n == 1
        omega(n) = 0;
        I_a(n) = 0;
    end
        omega(n+1) = (delta_t/J)*(K*I_a(n) - beta*omega(n)) + omega(n);
        I_a(n+1) = (delta_t/L) *(Va(n) - K*omega(n) - R*I_a(n)) + I_a(n);
end

averageval = mean(omega(end - 10000:end));
p2pval = peak2peak(omega(end - 10000:end));

subplot(2, 1, 1);
plot(t, Va);
xlabel('Time (secs)');
ylabel('Voltage (Volts)');
title('Voltage Source');

subplot(2, 1, 2);
plot(t, omega);
xlabel('Time (secs)');
ylabel('Speed (rads/s)');
title('Dependent Voltage Source');
