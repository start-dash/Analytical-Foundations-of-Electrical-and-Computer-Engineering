%% Project 2 - Problem 4

load('touch.mat');


%% part a - c
Fs = 8192; %% sampled at 8192 Hz
t = [0:1:999]; %% time interval
sample = 100; %% sample
deltat = 1/Fs; 
n_discretized = t*deltat; %% discretized time

d1 = sin(2*pi*697*n_discretized) + sin(2*pi*1209*n_discretized); %% 1
d2 = sin(2*pi*697*n_discretized) + sin(2*pi*1336*n_discretized); %% 2
d3 = sin(2*pi*697*n_discretized) + sin(2*pi*1477*n_discretized); %% 3
d4 = sin(2*pi*770*n_discretized) + sin(2*pi*1209*n_discretized); %% 4
d5 = sin(2*pi*770*n_discretized) + sin(2*pi*1336*n_discretized); %% 5
d6 = sin(2*pi*770*n_discretized) + sin(2*pi*1477*n_discretized); %% 6
d7 = sin(2*pi*852*n_discretized) + sin(2*pi*1209*n_discretized); %% 7
d8 = sin(2*pi*852*n_discretized) + sin(2*pi*1336*n_discretized); %% 8
d9 = sin(2*pi*852*n_discretized) + sin(2*pi*1477*n_discretized); %% 9
d0 = sin(2*pi*947*n_discretized) + sin(2*pi*1336*n_discretized); %% 0

space = zeros(size(n_discretized)); %% vector

x = [d7 space d0 space d4 space d4 space d2 space d1 space d3 space d4 space d4 space d6];
sound(x, Fs);



%% part d
N_p = 2048;
X = zeros(1, N_p);

d2p = [d2 zeros(1,2048 - length(d2))]; %% 1000 -> 2048 padding
d5p = [d5 zeros(1, 2048 - length(d5))]; %% 1000 -> 2048 padding

D2P = zeros(size(d2p));
D5P = zeros(size(d5p));

for k = 1:1:N_p  %% fourier transform
    for n = 1:1:N_p
            D2P(k) = D2P(k) + d2p(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            D5P(k) = D5P(k) + d5p(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
    end
end

k = [0:1:N_p - 1];
omega = k * 2 * pi/N_p;
w = k * Fs/N_p; %% convert radians a sec to Hz





%% part e
X1 = x1(1:1000);
X2 = x1(1101:2100);
X3 = x1(2201:3200);
X4 = x1(3301:4300);
X5 = x1(4401:5400);
X6 = x1(5501:6500);
X7 = x1(6601:7600);

DX1P = [X1 zeros(1,2048 - length(X1))]; %% 1000 -> 2048 padding
DX2P = [X2 zeros(1, 2048 - length(X2))]; %% 1000 -> 2048 padding
DX3P = [X3 zeros(1, 2048 - length(X3))]; %% 1000 -> 2048 padding
DX4P = [X4 zeros(1, 2048 - length(X4))]; %% 1000 -> 2048 padding
DX5P = [X5 zeros(1, 2048 - length(X5))]; %% 1000 -> 2048 padding
DX6P = [X6 zeros(1, 2048 - length(X6))]; %% 1000 -> 2048 padding
DX7P = [X7 zeros(1, 2048 - length(X7))]; %% 1000 -> 2048 padding

dx1p = zeros(size(DX1P));
dx2p = zeros(size(DX2P));
dx3p = zeros(size(DX3P));
dx4p = zeros(size(DX4P));
dx5p = zeros(size(DX5P));
dx6p = zeros(size(DX6P));
dx7p = zeros(size(DX7P));


for k = 1:1:N_p  %% fourier transform
    for n = 1:1:N_p
            dx1p(k) = dx1p(k) + DX1P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx2p(k) = dx2p(k) + DX2P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx3p(k) = dx3p(k) + DX3P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx4p(k) = dx4p(k) + DX4P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx5p(k) = dx5p(k) + DX5P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx6p(k) = dx6p(k) + DX6P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
            dx7p(k) = dx7p(k) + DX7P(n) * exp(-j * (2*pi)/N_p * (k-1) * (n-1));
    end
end


subplot(4, 2, 1);
plot(w, abs(dx1p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 1st Digit'); %% makes the number 4


subplot(4, 2, 2);
plot(w, abs(dx2p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 2nd Digit'); %% makes the number 9


subplot(4, 2, 3);
plot(w, abs(dx3p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 3rd Digit'); %% makes the number 1


subplot(4, 2, 4);
plot(w, abs(dx4p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 4th Digit'); %% makes the number 5


subplot(4, 2, 5);
plot(w, abs(dx5p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 5th Digit'); %% makes the number 8


subplot(4, 2, 6);
plot(w, abs(dx6p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 6th Digit'); %% makes the number 7


subplot(4, 2, 7);
plot(w, abs(dx7p));
xlim([0, 2048]);
xlabel('Omega');
ylabel('Magnitude');
title('Discrete Fourier Transform of the 7th Digit'); %% makes the number 7
