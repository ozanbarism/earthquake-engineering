clear 
clc

%An example linear elastic SDOF system is given with k = 4 kN/m and m = 15 kg.
%The system is initially at rest when it is excited by a harmonic force 
% F(t)=60sin12t (in Newtons).

%The first part of this code plots the curves for transient response, steady
%state response and total response separately for each damping ratio 

%The second part of this code plots response amplitude vs frequency ratio graphs 
%for each damping ratio considering the steady state response only and assuming
%that the excitation frequency 𝜔̅ is variable.

t=0:0.01:4; %s
%structural parameters
k=4000; %N/m
m=15; %kg
wn=sqrt(k/m);
%excitation parameters
F0=60; %Newtons
w_bar=12;
%damping ratios
beta=[0, 0.1, 0.2];


for j=1:length(beta)

wd(j)=wn*sqrt(1-beta(j)^2);
eps=w_bar/wn;
G1(j)=(F0/k)*((1-eps^2)/((1-eps^2)^2+(2*beta(j)*eps)^2));
G2(j)=(F0/k)*((-2*beta(j)*eps)/((1-eps^2)^2+(2*beta(j)*eps)^2));
%initial conditions assumed as displacement and velocity is zero at t=0;
%with this consideration, general solution is solved to find A1 and A2
A2(j)=(F0/k)*((2*beta(j)*eps)/((1-eps^2)^2+(2*beta(j)*eps)^2));
A1(j)=(((-F0*(1-eps^2)*w_bar)/(k*((1-eps^2)^2+(2*beta(j)*eps)^2)))+beta(j)*wn*A2(j))/wd(j);


    
    for i=1:length(t);
    uh(i,j)=exp(-beta(j)*wn*t(i))*(A1(j)*sin(wd(j)*t(i))+A2(j)*cos(wd(j)*t(i)));
    up(i,j)=G1(j)*sin(w_bar*t(i))+G2(j)*cos(w_bar*t(i));
    ug(i,j)=uh(i,j)+up(i,j);
    end
    


teta=atan(2*beta(j)*eps/(1-eps^2));
e=0:0.1:2;
w_bar_b=wn*e;
for i=1:length(e);
 ro(i,j)=(F0/k)/(sqrt((1-e(i).^2).^2+(2*beta(j)*e(i))^2));
end

end


figure();
plot(t,ug);
xlabel('t(s)');
ylabel('ug(m)');
legend('Damping Ratio=0%','Damping Ratio=10%','Damping Ratio=20%')

figure();
plot(t,uh);
xlabel('t(s)');
ylabel('uh(m))');
legend('Damping Ratio=0%','Damping Ratio=10%','Damping Ratio=20%')

figure();
plot(t,up);
xlabel('t(s)');
ylabel('up(m)');
legend('Damping Ratio=0%','Damping Ratio=10%','Damping Ratio=20%')

figure();
plot(e,ro);
xlabel('e');
ylabel('ro');
legend('Damping Ratio=0%','Damping Ratio=10%','Damping Ratio=20%')