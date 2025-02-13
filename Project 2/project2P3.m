%% Project 2 - Problem 3
clear all

load('problem3.mat');

%% part a
f = 60;
T = 1/f;
omega = 2*pi/T;

fs = 60000;
Ts = 1/fs;
t = [0:length(x)-1]*Ts;


N =[-13:1:13];
alpha_1 = 0;

for n = 1:1:length(t)
    alpha_1 = alpha_1 + x(n)*exp(-j*1*omega*t(n))*Ts;
end

alpha_1 = alpha_1/T;

figure(1);
plot(t, x);
title('Raw Input Signal');
xlabel('Time (sec)');
ylabel('Volts');

