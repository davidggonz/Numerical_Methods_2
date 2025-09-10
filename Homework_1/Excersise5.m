% Solve elliptic PDE: u_xx + u_yy = 4
clearvars;
% Domain: 0 < x < 1, 0 < y < 2
h = 1/2;
k = 1/2;

% Grid points
x = 0:h:1;
y = 0:k:2;

% Initialize solution with boundary conditions
U = zeros(length(y), length(x));
U(1,:) = x.^2;          % u(x,0) = x^2
U(end,:) = (x-2).^2;    % u(x,2) = (x-2)^2
U(:,1) = y.^2;          % u(0,y) = y^2
U(:,end) = (y-1).^2;    % u(1,y) = (y-1)^2

% Finite difference coefficients
alpha = 1/h^2;
beta = 1/k^2;
gamma = -2*(alpha + beta);

% Matrix A (3x3 system for interior points)
A = gamma * eye(3);

% Vector b (source term + boundary contributions)
b = 4 * ones(3,1);
b(1) = b(1) - beta*U(1,2) - alpha*U(2,1) - alpha*U(2,3);
b(2) = b(2) - beta*U(2,2) - alpha*U(3,1) - alpha*U(3,3);
b(3) = b(3) - beta*U(3,2) - alpha*U(4,1) - alpha*U(4,3);

% Solve
u_interior = A \ b;

% Fill interior points
U(2,2) = u_interior(1);
U(3,2) = u_interior(2);
U(4,2) = u_interior(3);

% Actual solution
U_actual = zeros(length(y), length(x));
for i = 1:length(y)
    for j = 1:length(x)
        U_actual(i,j) = (x(j) - y(i))^2;
    end
end

% Plot numerical solution vs x for each y
figure(1);
for i = 1:length(y)
    plot(x, U(i,:), '-o', 'DisplayName', sprintf('y = %.1f', y(i)));
end
title('Numerical Solution');
xlabel('x');
ylabel('u(x,y)');
legend;
grid on;

% Plot actual solution vs x for each y
figure(2);
hold on;
for i = 1:length(y)
    plot(x, U_actual(i,:), '-s', 'DisplayName', sprintf('y = %.1f', y(i)));
end
title('Actual Solution');
xlabel('x');
ylabel('u(x,y)');
legend;
grid on;
