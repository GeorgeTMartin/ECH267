%ECH 267 Homework 3
%George Martin
clc;clear;close;
global c count
t0 = 0;                     %initial time (almost always 0)
tf = 1;                     %final time (2s)
Nt = (tf-t0)*1000;          %number of total data points (1000 per second)
tspan = linspace(t0,tf,Nt);         %create linearly spaced data points from t0 to tf of size Nt
x0=[linspace(-5,5,11)];
c = 1;
for count = 1:2
    for i = 1:length(x0)
        [t,x] = ode45(@linfunc,tspan,x0(i)); 
        %extract outputs
        for j = 1:length(t)
            [dx(j,:)] = linfunc(t(i),x(j,:));
        end
        hold on
        subplot(2,1,count)
        plot(t,x)
        title(['X(t) - Method : ' num2str(count)])
        xlabel('t')
        ylabel('x')
    end
end
function [dx, dx2] = linfunc(t,x)
global c count
u(count) = -x*(1+c);
u(count) = - x - x * sqrt(1+x^4);
dx = x^3 + x^2 * u(count);
end
