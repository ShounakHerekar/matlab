%Code to find the specifications of a chebyshev type 2 
clear all;close all;clc;
Ap=2;
As=60;

wp=3000;
ws=7000;

e=1/sqrt((10^(As/10))-1);   %ripple factor

R=e^2/(1+e^2);          %ripple height

N=(acosh(sqrt((10^(As/10))-1)/sqrt(10^(Ap/10)-1)))/acosh(ws/wp);  %order of filter
N=ceil(N);  

wc=wp*cosh((1/N)*acosh(1/e))       %cutoff
[n,Ws] = cheb2ord(wp,ws,Ap,As,'s')

[b,a]=cheby2(N,As,ws,'s');
H=tf(b,a);      %transfer function
[h,w]=freqs(b,a);

figure(1)
plot(w,20*log10(abs(h)));       %magnitude response
xlabel('Frequency');
ylabel('Magnitude');
title('Magnitude response of Chebyshev type 2 filter')
grid on;
hold on;

figure(2)
pzmap(H)        %pole zero plot
grid on;
hold on;
