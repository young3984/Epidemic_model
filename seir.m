clear; clc; close all;

t=100;      % for 100 days, h=1day
beta=2; gamma=0.001; a=0.01;

S=zeros(1,t); E=S; I=S; R=S;
N=10010;
S(1)=10000/N; E(1)=10/N; I(1)=0; R(1)=0;

for ik=2:100
    S(ik)=S(ik-1)-beta*S(ik-1)*I(ik-1);
    E(ik)=E(ik-1)+beta*S(ik-1)*I(ik-1)-a*E(ik-1);
    I(ik)=I(ik-1)+a*E(ik-1)-gamma*I(ik-1);
    R(ik)=gamma*I(ik-1);
end

plot(S,'linewidth',2)
hold on
plot(E,'linewidth',2)
hold on
plot(I,'linewidth',2)
hold on
plot(R,'linewidth',2)
hold on