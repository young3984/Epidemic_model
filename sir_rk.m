% SIR: Runge Kutta
clear; close all; clc;

t=120;      % day
re=0.001;   % recover rate per time step
dt=t*re;
Brn=3.6;    % Basic reproductive number

N=10010;
S=zeros(1,t); I=S; R=S;

S(1)=10000/N; I(1)=10/N; R(1)=0/N;

for ik=2:t
    
    sh=S(ik-1); ih=I(ik-1);
    k1=dt*[-Brn*sh*ih Brn*sh*ih-ih ih];
    
    sh=(S(ik-1)+k1(1)*dt/2); ih=(I(ik-1)+k1(2)*dt/2);
    k2=dt*[-Brn*sh*ih Brn*sh*ih-ih ih];
    
    sh=(S(ik-1)+k2(1)*dt/2); ih=(I(ik-1)+k2(2)*dt/2);
    k3=dt*[-Brn*sh*ih Brn*sh*ih-ih ih];
    
    sh=(S(ik-1)+k3(1)); ih=(I(ik-1)+k3(2));
    k4=dt*[-Brn*sh*ih Brn*sh*ih-ih ih];
    
    S(ik)=S(ik-1)+(k1(1)+k2(1)+k3(1)+k4(1))/6;
    I(ik)=I(ik-1)+(k1(2)+k2(2)+k3(2)+k4(2))/6;
    R(ik)=R(ik-1)+(k1(3)+k2(3)+k3(3)+k4(3))/6;
end    

plot(S,'LineWidth',2)
hold on
plot(I,'LineWidth',2)
hold on 
plot(R,'LineWidth',2)
title("N=10010, S=10000, E=5, I=5, R=0, Brn=3.6, 100days, recover rate=0.001")
legend("S: Susceptible","I: Infectious","R: Recovered")
saveas(gcf,"SIR_RK4.jpg")