clear
clc

%extract data from text file
a_g=importdata('KLM86X01.txt');

%constant average acceleration variation is used
beta=1/4;
gama=1/2;
dt=0.01;
m=1;
%boundary conditions for three damping cases
v_1(1,1)=0; %initial velocity is zero
d_1(1,1)=0; %initial displacement is zero
v_2(1,1)=0;%initial velocity is zero
d_2(1,1)=0; %initial displacement is zero
v_3(1,1)=0;%initial velocity is zero
d_3(1,1)=0; %initial displacement is zero
%damping=0.02
b=0.02;
j=1;
%for loop for computing newmark method for every T value.
%inital T is taken as 0.0001 so that first values will not be infinity.
for T=0.0001:0.001:4
    k_1=(2*pi/T)^2*m;
    c_1=2*b*m*(2*pi/T);
    a_1(1,1)=(1/m)*(-m*a_g(1,1)-c_1*v_1(1,1)-d_1(1,1)*k_1);
    k_tilda_1=k_1+m*(1/(beta*dt^2))+c_1*(gama/(beta*dt));
    %for loop for newmark method
    for i=1:2585
        A_1=m*(1/(beta*dt))+c_1*gama/beta;
        B_1=m*1/(2*beta)+c_1*dt*(gama/(2*beta)-1);
        d_Ftilda_1(i,1)=(-m*a_g(i+1,1)+m*a_g(i,1))+A_1*v_1(i,1)+B_1*a_1(i,1);
        
        d_disp1(i,1)=d_Ftilda_1(i,1)/k_tilda_1;
        d_velo1(i,1)=(gama/(beta*dt))*d_disp1(i,1)-(gama/beta)*v_1(i,1)+dt*(1-gama/(2*beta))*a_1(i,1);
        d_acc1(i,1)=(1/(beta*dt^2))*d_disp1(i,1)-(1/(beta*dt))*v_1(i,1)-1/(2*beta)*a_1(i,1);
        
        d_1(i+1,1)=d_disp1(i,1)+d_1(i,1);
        v_1(i+1,1)=d_velo1(i,1)+v_1(i,1);
        a_1(i+1,1)=d_acc1(i,1)+a_1(i,1);
    end
    %extracting spectural values for each T
    Sd_1(j,1)=max(abs(d_1));
    Sv_1(j,1)=max(abs(v_1));
    Sa_1(j,1)=max(abs(a_1));
    Psv_1(j,1)=(2*pi/T)*Sd_1(j,1);
    Psa_1(j,1)=(2*pi/T)^2*Sd_1(j,1);
    j=j+1;
end

%same procedure done for damping=0.05
b=0.05;
j=1;
for T=0.0001:0.001:4
    k_2=(2*pi/T)^2*m;
    c_2=2*b*m*(2*pi/T);
    a_2(1,1)=(1/m)*(-m*a_g(1,1)-c_2*v_2(1,1)-d_2(1,1)*k_2);
    k_tilda_2=k_2+m*(1/(beta*dt^2))+c_2*(gama/(beta*dt));
    
    for i=1:2585
        A_2=m*(1/(beta*dt))+c_2*gama/beta;
        B_2=m*1/(2*beta)+c_2*dt*(gama/(2*beta)-1);
        d_Ftilda_2(i,1)=(-m*a_g(i+1,1)+m*a_g(i,1))+A_2*v_2(i,1)+B_2*a_2(i,1);
        
        d_disp2(i,1)=d_Ftilda_2(i,1)/k_tilda_2;
        d_velo2(i,1)=(gama/(beta*dt))*d_disp2(i,1)-(gama/beta)*v_2(i,1)+dt*(1-gama/(2*beta))*a_2(i,1);
        d_acc2(i,1)=(1/(beta*dt^2))*d_disp2(i,1)-(1/(beta*dt))*v_2(i,1)-1/(2*beta)*a_2(i,1);
        
        d_2(i+1,1)=d_disp2(i,1)+d_2(i,1);
        v_2(i+1,1)=d_velo2(i,1)+v_2(i,1);
        a_2(i+1,1)=d_acc2(i,1)+a_2(i,1);
    end
    Sd_2(j,1)=max(abs(d_2));
    Sv_2(j,1)=max(abs(v_2));
    Sa_2(j,1)=max(abs(a_2));
    Psv_2(j,1)=(2*pi/T)*Sd_2(j,1);
    Psa_2(j,1)=(2*pi/T)^2*Sd_2(j,1);
    j=j+1;
end

%same procedure done for damping=0.10
b=0.10;
j=1;
for T=0.0001:0.001:4
    k_3=(2*pi/T)^2*m;
    c_3=2*b*m*(2*pi/T);
    a_3(1,1)=(1/m)*(-m*a_g(1,1)-c_3*v_3(1,1)-d_3(1,1)*k_3);
    k_tilda_3=k_3+m*(1/(beta*dt^2))+c_3*(gama/(beta*dt));
    
    for i=1:2585
        A_3=m*(1/(beta*dt))+c_3*gama/beta;
        B_3=m*1/(2*beta)+c_3*dt*(gama/(2*beta)-1);
        d_Ftilda_3(i,1)=(-m*a_g(i+1,1)+m*a_g(i,1))+A_3*v_3(i,1)+B_3*a_3(i,1);
        
        d_disp3(i,1)=d_Ftilda_3(i,1)/k_tilda_3;
        d_velo3(i,1)=(gama/(beta*dt))*d_disp3(i,1)-(gama/beta)*v_3(i,1)+dt*(1-gama/(2*beta))*a_3(i,1);
        d_acc3(i,1)=(1/(beta*dt^2))*d_disp3(i,1)-(1/(beta*dt))*v_3(i,1)-1/(2*beta)*a_3(i,1);
        
        d_3(i+1,1)=d_disp3(i,1)+d_3(i,1);
        v_3(i+1,1)=d_velo3(i,1)+v_3(i,1);
        a_3(i+1,1)=d_acc3(i,1)+a_3(i,1);
    end
    Sd_3(j,1)=max(abs(d_3));
    Sv_3(j,1)=max(abs(v_3));
    Sa_3(j,1)=max(abs(a_3));
    Psv_3(j,1)=(2*pi/T)*Sd_3(j,1);
    Psa_3(j,1)=(2*pi/T)^2*Sd_3(j,1);
    j=j+1;
end

T=0.0001:0.001:4;

%graphs are plotted for different damping values
figure();
plot(T,Sa_1);
hold on
plot(T,Sa_2,'m');
hold on
plot(T,Sa_3,'r');
legend('Damping Ratio=2%','Damping Ratio=5%','Damping Ratio=10%')
xlabel('T(s)');
ylabel('Sa(m/s^2)');

figure();
plot(T,Sv_1);
hold on
plot(T,Sv_2,'m');
hold on
plot(T,Sv_3,'r');
legend('Damping Ratio=2%','Damping Ratio=5%','Damping Ratio=10%')
xlabel('T(s)');
ylabel('Sv(m/s)');

figure();
plot(T,Sd_1);
hold on
plot(T,Sd_2,'m');
hold on
plot(T,Sd_3,'r');
legend('Damping Ratio=2%','Damping Ratio=5%','Damping Ratio=10%')
xlabel('T(s)');
ylabel('Sd(m)');

figure();
plot(T,Psa_1);
hold on
plot(T,Psa_2,'m');
hold on
plot(T,Psa_3,'r');
legend('Damping Ratio=2%','Damping Ratio=5%','Damping Ratio=10%')
xlabel('T(s)');
ylabel('Psa(m/s^2)');


figure();
plot(T,Psv_1);
hold on
plot(T,Psv_2,'m');
hold on
plot(T,Psv_3,'r');
legend('Damping Ratio=2%','Damping Ratio=5%','Damping Ratio=10%')
xlabel('T(s)');
ylabel('Psv(m/s)');

    
