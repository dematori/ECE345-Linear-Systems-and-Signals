%% Jeffrey Huang - Linear Systems and Signals Lab 1 - Problem 1
%%
clc; clear all; close all;

syms e;
t = linspace(-1,1,1001);
u = @(t) (t>=0);
%% Dirac delta function 1
P1_e = @(e) (1./e*(u(t+(1/2)*e)-u(t-(1/2)*e)));
figure;
plot(t, P1_e(.05), 'b', t, P1_e(.1),'r:');
axis([-1 1 0 25]);
title('rectangular pulse');
xlabel('t');
legend('\bf\epsilon = 0.05', '\bf\epsilon = 0.10');
%% Dirac delta function 2
P2_e = @(e) (1./(sqrt(2*pi*e))*exp((-t.^2)/(2*e)));
figure;
plot(t, P2_e(0.001), 'b', t , P2_e(0.002), 'r:');
axis([-1 1 0 15]);
title('Gaussian');
xlabel('t');
legend('\bf\epsilon = 0.001', '\bf\epsilon = 0.002');
%% Dirac delta function 3
P3_e = @(e) ((1/pi)*(e./(e.^2+t.^2)));
figure;
plot(t, P3_e(0.1), 'b', t , P3_e(0.05), 'r:');
axis([-1 1 0 10]);
title('Gaussian');
xlabel('t');
legend('\bf\epsilon = 0.1', '\bf\epsilon = 0.05');
%% Dirac delta function 4
P4_e = @(e) (sinc(t/e)./(pi*t));
figure;
plot(t, P4_e(.001), 'b', t , P4_e(.002), 'r:');
title('Sinusoidal');
xlabel('t');
legend('\bf\epsilon = 0.001', '\bf\epsilon = 0.002');







