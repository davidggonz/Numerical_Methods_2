% Part (a): y'' + y = 0, y(0)=1, y(pi/2)=1
clearvars;
a = 0;
b = pi/2;
alpha = 1;
beta = 1;
N = 1000;


p = @(x) 0;
q = @(x) -1;
r = @(x) 0;

[x_a, w_a] = finite_differences(a, b, alpha, beta, N, p, q, r);

% Actual solution and error
y_actual_a = @(x) cos(x) + (sqrt(2)-1)*sin(x);
actual_vals_a = arrayfun(y_actual_a, x_a);
error_a = abs(w_a - actual_vals_a);

% Plot error
figure;
subplot(2,2,1);
plot(x_a, error_a, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Error');
title('(a) y'''' + y = 0');
grid on;
max_error_a = max(error_a);
mean_error_a = mean(error_a);
text(0.3, 0.2, sprintf('Max error: %.2e\nMean error: %.2e', max_error_a, mean_error_a),...
     'Units', 'normalized', 'BackgroundColor', 'white');

% Part (b): y'' + 4y = cos x, y(0)=0, y(pi/4)=0
a = 0;
b = pi/4;
alpha = 0;
beta = 0;
N = 1000;

p = @(x) 0;
q = @(x) -4;
r = @(x) cos(x);

[x_b, w_b] = finite_differences(a, b, alpha, beta, N, p, q, r);

% Actual solution and error
y_actual_b = @(x) -1/3*cos(2*x) - sqrt(6)/5*sin(2*x) + 1/3*cos(x);
actual_vals_b = arrayfun(y_actual_b, x_b);
error_b = abs(w_b - actual_vals_b);

% Plot error
subplot(2,2,2);
plot(x_b, error_b, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('Error');
title('(b) y'''' + 4y = cos x');
grid on;
max_error_b = max(error_b);
mean_error_b = mean(error_b);
text(0.3, 0.2, sprintf('Max error: %.2e\nMean error: %.2e', max_error_b, mean_error_b),...
     'Units', 'normalized', 'BackgroundColor', 'white');

% Part (c): y'' = -4x⁻¹y' + 2x⁻²y - 2x⁻²ln x, y(1)=0.5, y(2)=ln(2)
a = 1;
b = 2;
alpha = 0.5;
beta = log(2);
N = 1000;

p = @(x) -4/x;
q = @(x) 2/(x^2);
r = @(x) -2*log(x)/(x^2);

[x_c, w_c] = finite_differences(a, b, alpha, beta, N, p, q, r);

% Actual solution and error
y_actual_c = @(x) 4/x - 2/(x^2) + log(x) - 1.5;
actual_vals_c = arrayfun(y_actual_c, x_c);
error_c = abs(w_c - actual_vals_c);

% Plot error
subplot(2,2,3);
plot(x_c, error_c, 'g-', 'LineWidth', 2);
xlabel('x');
ylabel('Error');
title('(c) y'''' = -4x^{-1}y'' + 2x^{-2}y - 2x^{-2}ln x');
grid on;
max_error_c = max(error_c);
mean_error_c = mean(error_c);
text(0.3, 0.2, sprintf('Max error: %.2e\nMean error: %.2e', max_error_c, mean_error_c),...
     'Units', 'normalized', 'BackgroundColor', 'white');

% Part (d): y'' = 2y' - y + xeˣ - x, y(0)=0, y(2)=-4
a = 0;
b = 2;
alpha = 0;
beta = -4;
N = 1000;

p = @(x) 2;
q = @(x) -1;
r = @(x) x.*exp(x) - x;

[x_d, w_d] = finite_differences(a, b, alpha, beta, N, p, q, r);

% Actual solution and error
y_actual_d = @(x) 1/6*x.^3.*exp(x) - 3/2*x.*exp(x) + 2*exp(x) - x - 2;
actual_vals_d = arrayfun(y_actual_d, x_d);
error_d = abs(w_d - actual_vals_d);

% Plot error
subplot(2,2,4);
plot(x_d, error_d, 'm-', 'LineWidth', 2);
xlabel('x');
ylabel('Error');
title('(d) y'''' = 2y'' - y + xe^x - x');
grid on;
max_error_d = max(error_d);
mean_error_d = mean(error_d);
text(0.3, 0.2, sprintf('Max error: %.2e\nMean error: %.2e', max_error_d, mean_error_d),...
     'Units', 'normalized', 'BackgroundColor', 'white');