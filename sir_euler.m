% SIR: Euler Method
clear; close all; clc;

t=100;      % day
re=0.001;   % recover rate per time step
dt=t*re;
Brn=3.6;    % Basic reproductive number

N=10010;
S=zeros(1,t); I=N; R=N;

S(1)=10000/N; I(1)=10/N; R(1)=0/N;

for ik=2:t
    S(ik)=S(ik-1)+dt*(-Brn*S(ik-1)*I(ik-1));
    I(ik)=I(ik-1)+dt*(Brn*S(ik-1)*I(ik-1)-I(ik-1));
    R(ik)=R(ik-1)+dt*I(ik-1);   
end    

plot(S,'LineWidth',2)S
hold on
plot(I,'LineWidth',2)
hold on 
plot(R,'LineWidth',2)
title("N=10010, S=10000, E=5, I=5, R=0, Brn=3.6, 100days, recover rate=0.001")
legend("S: Susceptible","I: Infectious","R: Recovered")
saveas(gcf,"SIR_Euler.jpg")