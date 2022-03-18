%ECH267 Final Project
%George Martin
clc;clear;close;
global K1 K2 K3 K4 K5 K6 K7 K8 T_amb T_wo tprev k Ierror MaxFlow Gain
rho_a = 1.225; %kg/m^3
cp_a = 1000; %J/[Kg*K]
M_tc = 25; %Thermal Mass Capacitance of Room
Q_w = .0005; %m^3/s [Ratio of Chiller Flow to Output Air]
c_pw = 4186; %J/[kg*K]
M_cc = .5; %kg
c_v = 718; %J/[kg*K]
V_t = 3; %m^3
rho_w = 1000; %kg/m^3
T_wo = 10; %C Temperature of Water Leaving Chiller
T_amb = 30; %C Ambient Temperature
MaxFlow = .5; %Maximum Fan Flow Rate
U_cc = 15; %Heat Transfer Coefficient of Chiller Coil
U_tc = 30; %Heat Transfer Coefficient of Room
U_ta = 10; %Heat Transfer Coefficient of Tank
A_ta = 2; %Surface Area of Tank
A_cc = 5; %Surface Area of Chiller Coil
Qgen = 150; %J Heat Produced in Room
Gain =.4;                                                       %Change Gain Value to Increase or Decrease Response

K1 = rho_a*cp_a/(M_tc*c_v);
K2 = U_tc*A_cc/(M_tc*c_v);
K3 = Qgen/(M_tc*c_v);
K4 = U_ta*A_ta/(V_t*rho_w*c_pw);
K5 = Q_w*rho_w*c_pw/(V_t*rho_w*c_pw);
K6 = rho_a*cp_a/(V_t*rho_w*c_pw);
K7 = U_cc*A_cc/(V_t*rho_w*c_pw);
K8 = Q_w*rho_w/(M_cc*c_v);

x10 = T_amb; %Room Temperature Start
x20 = 10; %Tank Temperature Start
x30 = 10; %Exchange Air Temperature Start

x0=[x10 x20 x30]';     %linear initial conditions
tspan = linspace(0,600,3000);

%Solver---------------------------------------------------
%LINEAR
figure 
hold on
for k = 1:3
tprev =0;
Ierror = 0; %Initialize Integral Error
[t,x] = ode45(@linfunc,tspan,x0);       %ode solver for linear case

for i = 1:length(t)
    [dx(i,:),extra(i,:)] = linfunc(t(i),x(i,:));
end
u = extra(:,1);
Perror = extra(:,2);
Ierror= extra(:,3);
Temp_room = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

subplot(3,1,1)
plot(t,Temp_room)
title('Room Temperature [C] v. Time [s]')
yline(25)
ylim([24 T_amb+5])
hold on
subplot(3,1,2)
plot(t,u)
title('Controlled Air Flow [m^3/s] v. Time [s]')
hold on
subplot(3,1,3)
plot(t,Perror)
title('Error [C] v. Time [s]')
hold on
SteadyStateMaximum = max(Temp_room(1000:3000))
end
legend('No Input','Linear','Nonlinear')

function [dx, extra] = linfunc(t,x)
global K1 K2 K3 K4 K5 K6 K7 K8 T_amb T_wo k Ierror t_prev MaxFlow Gain
x1 = x(1);
x2 = x(2);
x3 = x(3);
Perror = x1-25;
dt = t-t_prev;
%Ierror = Ierror + Perror*dt;
switch k
    case 1 %ZERO INPUT
        u = 0;
    case 2 %CLASSICAL PID CONTROLLER
        if x1 >25
        u = Gain*(Perror);
            if u > MaxFlow
                u=MaxFlow;
            end
        else
            u=0;
        end
    case 3
        %-.4*(K1*(x3-x1) + K2*(T_amb-x1)+K3)
         u = ((-Gain*Perror)-(K2*(T_amb-x1)+K3))/(K1*(x3-x1));
            if u > MaxFlow
             u = MaxFlow;
            end
        if u<0
            u=0;
        end
out1=K2*(T_amb-x1)+K3;
out2=K1*(x3-x1);
end
dx1 = K1*u*(x3-x1) + K2*(T_amb-x1)+K3; %T Room
dx2 = K4*(T_amb-x1)+K5*(T_wo-x2); %T Air Out
dx3 = K6*u*(x1-x3)+K7*(T_amb-x3+x2)+K8*u*(x2-T_wo); %T water tank
t_prev = t;
dx = [dx1 dx2 dx3]';
extra = [u Perror Ierror];

end
