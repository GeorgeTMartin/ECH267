% ECH 267 - Homework 1
% George Martin - 7 Feb 2022
clear all; close all; clc
% Problem 3 ----------------------------------------------------%
figure(1)
k = 5;a = .15;F = 2;
f = @(t,X) [X(2);-k*(1-a^2*X(1)^2)*X(1)-sign(X(1))*1.5+F];
[x,y] = meshgrid(-2:.1:2,-4:.1:4);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis tight equal;
hold on
for y20 = -1:5
    [ts,ys] = ode45(f,[0,10],[0;y20]);
    plot(ys(:,1),ys(:,2))
end
hold off

% Problem 6 ----------------------------------------------------%
figure(2)
for j = 1:1:5
if j==1
 A=[0 1;-2 -3]; Z=[-1 0;0 -2];
end
if j==2
 A=[0 -1;1 2]; Z=[1 0;0 1];
end
if j==3
 A=[1 1;0 -1]; Z=[1 0;0 -1];
end
if j==4
 A=[1 5;-1 -1]; Z=[0 2;-2 0];
end
if j==5
 A=[2 -1;2 0]; Z=[1 1;-1 1];
end
f = @(t,X) [A(1,1)*X(1)+A(1,2)*X(2);A(2,1)*X(1)+A(2,2)*X(2)];
f2 = @(t,X) [Z(1,1)*X(1)+Z(1,2)*X(2);Z(2,1)*X(1)+Z(2,2)*X(2)];
[x,y] = meshgrid(-2:.1:2,-2:.1:2);
u = zeros(size(x));
v = zeros(size(x));
u2 = zeros(size(x));
v2 = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    Zprime = f2(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
    u2(i) = Zprime(1);
    v2(i) = Zprime(2);
end
subplot(5,2,2*j-1)
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis tight equal;
hold on
for y20 = -1:.5:1
    for x20 = -1:.5:1
    [ts,ys] = ode45(f,[0,50],[x20;y20]);
    subplot(5,2,2*j-1)
    plot(ys(:,1),ys(:,2))
    axis([-2 2 -2 2]);
    end
end
hold off
subplot(5,2,2*j)
quiver(x,y,u2,v2,'r'); figure(gcf)
xlabel('X_1')
ylabel('X_2')
axis([-2 2 -2 2])
hold on
for y20 = -1:.5:1
    for x20 = -1:.5:1
    [ts,ys] = ode45(f2,[0,50],[x20;y20]);
    subplot(5,2,2*j)
    plot(ys(:,1),ys(:,2))
    end
end
hold off
end

% Problem 7 ------------------------------------------------------%
figure(3)
f = @(t,X) [-X(1)-X(2)/(log(sqrt(X(1)^2+X(2)^2)));-X(2)+X(1)/(log(sqrt(X(1)^2+X(2)^2)))];
[x,y] = meshgrid(-1:.1:1,-1:.1:1);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end
u(isinf(u)|isnan(u)) = 0;v(isinf(v)|isnan(v)) = 0;
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis tight equal;
hold on
for y20 = -.5:.25:.5
    for x20 = -.5:.25:.5
    [ts,ys] = ode45(f,[0,10],[x20;y20]);
    plot(ys(:,1),ys(:,2))
    end
end
hold off

% Problem 8 ------------------------------------------------------%
figure(4)
for j = 1:1:4
if j==1
  f = @(t,X) [X(2); X(1)-2*atan(X(1)+X(2))];
end
if j==2
  f = @(t,X) [X(2); -X(1)+X(2)*(1-3*X(1)^2-2*X(2)^2)]; 
end
if j==3
  f = @(t,X) [X(1)-X(1)*X(2); 2*X(1)^2-X(2)];
end
if j==4
  f = @(t,X) [X(1)+X(2)-X(1)*(abs(X(1)+abs(X(2)))); -2*X(1)+X(2)-X(2)*(abs(X(1)+abs(X(2))))];
end
[x,y] = meshgrid(-2:.1:2,-2:.1:2);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end
subplot(2,2,j)
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis tight equal;
hold on
for y20 = -1:.5:1
    for x20 = -1:.5:1
    [ts,ys] = ode45(f,[0,50],[x20;y20]);
    subplot(2,2,j)
    plot(ys(:,1),ys(:,2))
    end
end
hold off
end

% Problem 9 ------------------------------------------------
figure(5)
for j = 1:1:4
if j==1
  f = @(t,X) [-X(2); X(1)-X(2)*(1-X(1)^2+.1*X(1)^4)];
end
if j==2
  f = @(t,X) [X(2); X(1)+X(2)-3*atan(X(1)+X(2))]; 
end
if j==3
  f = @(t,X) [X(2); -(.5*X(1)+X(1)^3)];
end
if j==4
 old = sympref('HeavisideAtOrigin',1);
  f = @(t,X) [X(2); -X(2)-(heaviside(X(1)-X(2)-1)*(2*(X(1)-X(2))-.5*sign((X(1)-X(2))))+(heaviside(-(X(1)-X(2)-1))*(X(1)-X(2)^3+.5*(X(1)-X(2)))))];
end
[x,y] = meshgrid(-1.5:.1:1.5,-1.5:.1:1.5);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    Xprime = f(t,[x(i); y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end
subplot(2,2,j)
quiver(x,y,u,v,'r');
xlabel('X_1')
ylabel('X_2')
axis ([-1.5 1.5 -1.5 1.5]);
hold on
for y20 = -.5:.25:.5
    for x20 = -.5:.25:.5
    [ts,ys] = ode45(f,[0,10],[x20;y20]);
    subplot(2,2,j)
    plot(ys(:,1),ys(:,2))
    end
end
hold off
end

