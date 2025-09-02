%Matrix

A = [ [5, -2, 8, -7]; ...
      [12, 15, -14, -20];... 
      [8, 1, -3, -4];...
      [11, 18, 21, -3] ];
%Column vector transpouse
B =  [5, -2, 8, -7]';

%Inverse of a matrix
A_inv1 = inv(A); 
A_inv2 = A^(-1);

%Solutions:
x1 = A_inv1*B; % Gauss Jordan
x2 = inv(A) *B; % Gauss Jordan
x3 = A^(-1)*B; % Gauss Jordan
x4 = A\B; % Gauss-Jordan faster

figure(1)
spy(A_inv1)


I = A\ A; %Identity
figure(2)
spy(I)

%This depends entirely on the precision
N = 1e4;
b = rand(N);
c = rand(N,1);

%Solucions 
disp(' y1 = b^(-1)*c;');
tic;
y1 = b^(-1)*c;
toc;

disp(' y2 = inv(b)*c;');

tic;
y2 = inv(b)*c;
toc;

disp(' y3 = b\c;');

tic;
y3 = b\c;
toc;

% This shows that the last one is the fastest
%Determinant
disp('D = det(b); ');
tic;
D = det(b);
toc;

delta1 = y3-y2;
delta2 = y3-y1;
figure(1);
plot(delta1, '.r');
hold on;
plot(delta2, '.b');
hold off;