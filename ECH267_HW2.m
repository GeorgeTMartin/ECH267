% ECH 267 - Homework 2
% George Martin - 14 Feb 2022
clc;clear;close;
%PROBLEM 6B 
t0 = 0;                     %initial time (almost always 0)
tf = 5;                     %final time (2s)
Nt = (tf-t0)*1000;          %number of total data points (1000 per second)
tspan = linspace(t0,tf,Nt);         %create linearly spaced data points from t0 to tf of size Nt
x0=1;
[t,x] = ode45(@linfunc,tspan,x0); 
%extract outputs
for i = 1:length(t)
    [dx(i,:)] = linfunc(t(i),x(i,:));
    V(i)=x(i)^2;
end
subplot(2,1,1)
plot(t,x)
title('X(t)')
xlabel('t')
ylabel('x')
subplot(2,1,2)
plot(t,V)
title('V(x(t))')
xlabel('t')
ylabel('V')

%PROBLEM 10
f = @(t,X) [X(2); -X(1)+(1/3)*X(1)^3-X(2)];
[x,y] = meshgrid(-2:.1:2,-2:.1:2);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end
figure
subplot(2,1,1)
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis ([-2 2 -2 2]);
hold on
for y20 = -2:.5:2
    for x20 = -2:.5:2
    [ts,ys] = ode45(f,[0,4],[x20;y20]);
    plot(ys(:,1),ys(:,2))
    end
end

N=0:.125:2;
V=.75.*x.^2-(1/12).*x.^4+.5.*x.*y+.5.*y.^2;
[C,h]=contour(x,y,V,[N(10) N(10)]);
clabel(C,h)
h.LineWidth = 3;
hold off

subplot(2,1,2)
[C,h]=contour(x,y,V,N);
clabel(C,h)

%FUNCTION FOR PROBLEM 6B
function [dx] = linfunc(t,x)
dx=-x^3;
end