%Excersises 2-4
function nonlinear_finite_difference()
    % Part (a): y'' = y^3 - yy', 1 ≤ x ≤ 2, y(1)=1/2, y(2)=1/3
    disp('Part (a):');
    a = 1; b = 2; alpha = 1/2; beta = 1/3; h = 0.1;
    x = a + h:h:b - h;
    
    w0 = alpha + (beta - alpha)/(b - a) * (x - a);
    
    % Solve using Newton-Raphson 
    %Function, Jacobian, ω0, h, α, β tolerance, N
    [w_a, iter_a] = newton_raphson_nonlinear(@F_a, @J_a, w0, h, alpha, beta, 1e-8, 1000);
    
    y_actual_a = @(x) 1./(x + 1);
    x_full = [a, x, b];
    actual_vals_a = arrayfun(y_actual_a, x_full);
    numerical_vals_a = [alpha; w_a(:); beta]; % Convert to column vector
    error_a = abs(numerical_vals_a - actual_vals_a(:));
    
    figure;
    subplot(3,2,1);
    plot(x_full, numerical_vals_a, 'b-o', 'LineWidth', 2, 'DisplayName', 'Numerical');
    hold on;
    xx = linspace(a, b, 100);
    plot(xx, y_actual_a(xx), 'r--', 'LineWidth', 2, 'DisplayName', 'Actual');
    title('(a) y'''' = y^3 - yy''');
    xlabel('x'); ylabel('y(x)'); legend; grid on;
    
    subplot(3,2,2);
    plot(x_full, error_a, 'm-', 'LineWidth', 2);
    title('Error'); xlabel('x'); ylabel('|Error|'); grid on;
    
    % Part (b): y'' = 2y^3 - 6y - 2x^3, 1 ≤ x ≤ 2, y(1)=2, y(2)=5/2
    disp('Part (b):');
    a = 1; b = 2; alpha = 2; beta = 2.5; h = 0.1;
    N = 1000;
    x = a + h:h:b - h;
    
    w0 = alpha + (beta - alpha)/(b - a) * (x - a);
    
    %Function, Jacobian, ω0, h, α, β tolerance, N
    [w_b, iter_b] = newton_raphson_nonlinear(@F_b, @J_b, w0, h, alpha, beta, 1e-8, 1000);
    
    y_actual_b = @(x) x + 1./x;
    x_full = [a, x, b];
    actual_vals_b = arrayfun(y_actual_b, x_full);
    numerical_vals_b = [alpha; w_b(:); beta];
    error_b = abs(numerical_vals_b - actual_vals_b(:));
    
    subplot(3,2,3);
    plot(x_full, numerical_vals_b, 'b-o', 'LineWidth', 2, 'DisplayName', 'Numerical');
    hold on;
    xx = linspace(a, b, 100);
    plot(xx, y_actual_b(xx), 'r--', 'LineWidth', 2, 'DisplayName', 'Actual');
    title('(b) y'''' = 2y^3 - 6y - 2x^3');
    xlabel('x'); ylabel('y(x)'); legend; grid on;
    
    subplot(3,2,4);
    plot(x_full, error_b, 'm-', 'LineWidth', 2);
    title('Error'); xlabel('x'); ylabel('|Error|'); grid on;
    
    % Part (d): y'' = (y')^2x^{-3} - 9y^2x^{-5} + 4x, 1 ≤ x ≤ 2, y(1)=0, y(2)=ln(256)
    disp('Part (d):');
    a = 1; b = 2; alpha = 0; beta = log(256); h = 0.05;
    N = round((b - a)/h - 1);
    x = a + h:h:b - h;
    
    w0 = alpha + (beta - alpha)/(b - a) * (x - a);

    %Function, Jacobian, ω0, h, α, β tolerance, N
    [w_d, iter_d] = newton_raphson_nonlinear(@F_d, @J_d, w0, h, alpha, beta, 1e-8, 1000);
    
    y_actual_d = @(x) x.^3 .* log(x);
    x_full = [a, x, b];
    actual_vals_d = arrayfun(y_actual_d, x_full);
    numerical_vals_d = [alpha; w_d(:); beta];
    error_d = abs(numerical_vals_d - actual_vals_d(:));
    
    subplot(3,2,5);
    plot(x_full, numerical_vals_d, 'b-o', 'LineWidth', 2, 'DisplayName', 'Numerical');
    hold on;
    xx = linspace(a, b, 100);
    plot(xx, y_actual_d(xx), 'r--', 'LineWidth', 2, 'DisplayName', 'Actual');
    title('(d) y'''' = (y'')^2x^{-3} - 9y^2x^{-5} + 4x');
    xlabel('x'); ylabel('y(x)'); legend; grid on;
    
    subplot(3,2,6);
    plot(x_full, error_d, 'm-', 'LineWidth', 2);
    title('Error'); xlabel('x'); ylabel('|Error|'); grid on;
end



% Functions
function F = F_a(w, h, alpha, beta)
    N = length(w);
    F = zeros(N, 1);
    
    for i = 1:N
        if i == 1
            y_prev = alpha;
            y = w(1);
            y_next = w(2);
        elseif i == N
            y_prev = w(N-1);
            y = w(N);
            y_next = beta;
        else
            y_prev = w(i-1);
            y = w(i);
            y_next = w(i+1);
        end
        
        % Finite difference approximations
        y_prime = (y_next - y_prev) / (2*h);
        y_double_prime = (y_prev - 2*y + y_next) / (h^2);
        
        % Nonlinear equation: y'' - y^3 + y*y' = 0
        F(i) = y_double_prime - y^3 + y*y_prime;
    end
end

function J = J_a(w, h, alpha, beta)
    N = length(w);
    J = zeros(N, N);
    
    for i = 1:N
        if i == 1
            y_prev = alpha;
            y = w(1);
            y_next = w(2);
        elseif i == N
            y_prev = w(N-1);
            y = w(N);
            y_next = beta;
        else
            y_prev = w(i-1);
            y = w(i);
            y_next = w(i+1);
        end
        
        % Derivatives of the nonlinear equation
        dF_dy_prev = 1/h^2 - y/(2*h);
        dF_dy = -2/h^2 - 3*y^2 + (y_next - y_prev)/(2*h);
        dF_dy_next = 1/h^2 + y/(2*h);
        
        if i > 1
            J(i, i-1) = dF_dy_prev;
        end
        J(i, i) = dF_dy;
        if i < N
            J(i, i+1) = dF_dy_next;
        end
    end
end

function F = F_b(w, h, alpha, beta)
    N = length(w);
    F = zeros(N, 1);
    x = 1 + h:h:2 - h;  % x values for interior points
    
    for i = 1:N
        if i == 1
            y_prev = alpha;
            y = w(1);
            y_next = w(2);
        elseif i == N
            y_prev = w(N-1);
            y = w(N);
            y_next = beta;
        else
            y_prev = w(i-1);
            y = w(i);
            y_next = w(i+1);
        end
        
        y_double_prime = (y_prev - 2*y + y_next) / (h^2);
        F(i) = y_double_prime - 2*y^3 + 6*y + 2*x(i)^3;
    end
end

function J = J_b(w, h, alpha, beta)
    N = length(w);
    J = zeros(N, N);
    x = 1 + h:h:2 - h;
    
    for i = 1:N
        if i == 1
            y = w(1);
        elseif i == N
            y = w(N);
        else
            y = w(i);
        end
        
        dF_dy_prev = 1/h^2;
        dF_dy = -2/h^2 - 6*y^2 + 6;
        dF_dy_next = 1/h^2;
        
        if i > 1
            J(i, i-1) = dF_dy_prev;
        end
        J(i, i) = dF_dy;
        if i < N
            J(i, i+1) = dF_dy_next;
        end
    end
end

function F = F_d(w, h, alpha, beta)
    N = length(w);
    F = zeros(N, 1);
    x = 1 + h:h:2 - h;  % x values for interior points
    
    for i = 1:N
        if i == 1
            y_prev = alpha;
            y = w(1);
            y_next = w(2);
        elseif i == N
            y_prev = w(N-1);
            y = w(N);
            y_next = beta;
        else
            y_prev = w(i-1);
            y = w(i);
            y_next = w(i+1);
        end
        
        y_prime = (y_next - y_prev) / (2*h);
        y_double_prime = (y_prev - 2*y + y_next) / (h^2);
        
        F(i) = y_double_prime - (y_prime)^2/(x(i)^3) + 9*y^2/(x(i)^5) - 4*x(i);
    end
end

function J = J_d(w, h, alpha, beta)
    N = length(w);
    J = zeros(N, N);
    x = 1 + h:h:2 - h;
    
    for i = 1:N
        if i == 1
            y_prev = alpha;
            y = w(1);
            y_next = w(2);
        elseif i == N
            y_prev = w(N-1);
            y = w(N);
            y_next = beta;
        else
            y_prev = w(i-1);
            y = w(i);
            y_next = w(i+1);
        end
        
        y_prime = (y_next - y_prev) / (2*h);
        
        dF_dy_prev = 1/h^2 + (y_prime)/(x(i)^3 * h);
        dF_dy = -2/h^2 + 18*y/(x(i)^5);
        dF_dy_next = 1/h^2 - (y_prime)/(x(i)^3 * h);
        
        if i > 1
            J(i, i-1) = dF_dy_prev;
        end
        J(i, i) = dF_dy;
        if i < N
            J(i, i+1) = dF_dy_next;
        end
    end
end