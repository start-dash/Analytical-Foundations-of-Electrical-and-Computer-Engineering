%% Project 2 - Problem 2

%% part a
T = 10^-4; %% period, a.a
f = 1/T;
w = (2*pi*f); 
delta_t = T/10000;
t = [-2*T:delta_t:2*T]; %% time step, a.c
x = -sawtooth(2*pi/T*t, 0.5); %%% sawtooth command & phase, a.b & a.d

% plot(t, x);
% figure(1);
% title('Sawtooth Waveform');
% xlabel('Time: sec');
% grid on;



%% determine the corresponding signals, part a
alpha_1 = 0;
alpha_2 = 0;
alpha_3 = 0;
for n = 1:1:10000
     alpha_1 = alpha_1 + x(n)*exp(-j*1*w*t(n))*delta_t; %% create the waveform x
     alpha_2 = alpha_2 + x(n)*exp(-j*2*w*t(n))*delta_t;
     alpha_3 = alpha_3 + x(n)*exp(-j*3*w*t(n))*delta_t;
end

alpha_1 = alpha_1/T;
alpha_2 = alpha_2/T;
alpha_3 = alpha_3/T;
abs(alpha_1);
abs(alpha_2);
abs(alpha_3);




%% part b
h1w = @(wc,w) 1/(1+((j*w)/wc)); %% h1w takes in wc and w and outputs H1 @ omega
h2w = @(wc, w) (wc)^2./( (s).^2 + (s .* wc .* ( 2/sqrt(2) ) ) + (wc)^2 ); %% h2w takes in wc and w and outputs H2 @ omega
h4w = @(wc, w) (wc)^4./( ((s).^2 + (s .* wc .* 0.7654) + (wc)^2) .*  ( (s).^2 + (s.*wc*1.8478) + (wc)^2 )); %% h4w takes in wc and w and outputs H4 @ omega

wc = 2*pi*100; %%TO DO: play around with omega c and why the first 2 won't work and why 3rd one will.

w0 = (w);
w2pi60 = (2*pi*60);

hw0 = h1w(wc,w0);
hw2pi60=h1w(wc, w2pi60);

error = 0.01;


if abs( abs(hw0) - 0.01 ) / 0.01 <= error %% experimental value - expected value / expected value = percent error
    disp('passes const. 1');
else
    disp('fails const. 1');
end

if abs( angle(hw2pi60) ) <= 4 %% experimental value < 4
    disp('passes const. 2');
else
    disp('fails const. 2');
end

if abs( abs(hw2pi60) - 1 ) <= error
    disp('passes const. 3');
else
    disp('fails const. 3');
end

%% part c
omega = logspace(1,6,10000); %% go all the way to the frequency



s = j.*omega;
H1W = 1./(1+(s./wc)); %% H(1)W
H2W = (wc)^2./( (s).^2 + (s .* wc .* ( 2/sqrt(2) ) ) + (wc)^2 ); %% H(2)W
H4W = (wc)^4./( ((s).^2 + (s .* wc .* 0.7654) + (wc)^2) .*  ( (s).^2 + (s.*wc*1.8478) + (wc)^2 )); %% H(4)W


% %% H(1)W
% subplot(2, 1, 1);
% semilogx(omega, 20*log10(abs(H1W))); %% magnitude plot
% title('H1 Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');
% subplot(2, 1, 2);
% semilogx(omega, rad2deg(abs(H1W))); %% bode plot
% title('H1 Angle');
% xlabel('Frequency');
% ylabel('Angle');

% %% H(2)W
% subplot(2, 1, 1);
% semilogx(omega, 20*log10(abs(H2W))); %% magnitude plot
% title('H2 Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');
% subplot(2, 1, 2);
% semilogx(omega, rad2deg(abs(H2W))); %% angle plot
% title('H2 Angle');
% xlabel('Frequency');
% ylabel('Angle');

% %% H(4)W
% subplot(2, 1, 1);
% semilogx(omega, 20*log10(abs(H4W))); %% magnitude plot
% title('H4 Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');
% subplot(2, 1, 2);
% semilogx(omega, rad2deg(angle(H4W))) %% angle plot
% title('H4 Angle');
% xlabel('Frequency');
% ylabel('Angle');




%% part d, filter 1 and 2 angles immediately go down. you'll meet 2 constraints but not your 3rd
alpha_4 = 1; %% constructs the a_n values for the Vin signal

h3w = @ (w, wc) (wc)^4./( ((s).^2 + (s .* wc .* 0.7654) + (wc)^2) .*  ( (s).^2 + (s.*wc*1.8478) + (wc)^2 )); %% New filter

% construct alpha values for Vout
h1 = h3w(w, wc);
h2 = h3w(2*w, wc);
h3 = h3w(2*pi*60, wc);
h4 = h3w(3*w, wc);
% create the alpha values
v_a1 = (alpha_1 * h1);
v_a2 = (alpha_2 * h2);
v_a3 = (alpha_3 * h3);
v_a4 = (alpha_4 * h4);
% create V_out
V_out1 = 2 * abs(v_a1) * cos(w * t + angle(v_a1));
V_out2 = 2 * abs(v_a2) * cos(2* w * t + angle(v_a2));
V_out3 = 2 * abs(v_a3) * cos(3 * w * t + angle(v_a3));
V_out4 = 2 * abs(v_a4) * cos(2 * pi * 60 * t + angle(v_a4));

V_out = V_out1 + V_out2 + V_out3 + V_out4;


figure(6);
subplot(2, 1, 1);
plot(t, V_out);
title('Filtered Voltage Out');
xlabel('Time: sec');
ylabel('Voltage');
subplot(2, 1, 2);
plot(t, Vs);
title('Unfiltered Voltage Source');
xlabel('Time: sec');
ylabel('Voltage');
