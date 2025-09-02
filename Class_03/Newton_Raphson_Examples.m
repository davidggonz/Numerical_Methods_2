% Metodo Newton - Raphson
%Example: f(x)  = (x-2)(x+5)(x-7) y su derivada f'(x) 
%This show the zeros
clearvars;
x = (-6:0.001:10)';
%The function
f = @(x) (x-2).*(x+5).*(x-7);
%Its derivative
df = @(x) (x+5).*(x-7) + (x-2).*(x-7) + (x-2).*(x+5);
figure(1);
plot(x,f(x), '-r');
tol = 1e-10;
n = 100;
xn = 0;
for i = 1:n
    xnn = xn - f(xn)/df(xn);
    if abs(xnn-xn)<=tol
        break;
    end
    xn = xnn;
end
%Let's use it so we can find the roots for a non linear equations system
% Ejemplo 2: 4x^2 + y^2 + 2*x*y - y = 2 ====> 4x^2 + y^2 + 2*x*y - y -2= 0
%            2x^2 + 3*x*y + y^2 = 3     ====> 2x^2 + 3*x*y + y^2 -3 = 0

%Initial values
x = 0.1;
y = 0.1;

sn = [x,y]';
tol = 1e-10;
N = 100;
for i = 1:N
    J = [8*x+2*y, 2*y + 2*x-1;...
        4*x+3*y, 3*x + 2*y];%Jacobian
    F = [4*x^2+y^2+2*x*y-y-2;...
        2*x^2+3*x*y+y^2-3];%The system
    v = J\(-F); % 
    snn = sn + v;
    test = sqrt(dot(sn - snn, sn-snn));
    if test<= tol
        break;
    end
    sn = snn;
    x = sn(1);
    y = sn(2);
end