function [x, w] = finite_differences(a, b, alpha, beta, N, p, q, r,h)
    % Step 1: Initialize variables
    h = (b - a) / (N + 1);
    x = a + h;
    
    % Initialize arrays
    a_arr = zeros(N, 1);
    b_arr = zeros(N, 1);
    c_arr = zeros(N, 1);
    d_arr = zeros(N, 1);
    
    % Step 1: First equation (i=1)
    a_arr(1) = 2 + h^2 * q(x);
    b_arr(1) = -1 + (h/2) * p(x);
    d_arr(1) = -h^2 * r(x) + (1 + (h/2) * p(x)) * alpha;
    
    % Step 2: Interior equations (i=2 to N-1)
    for i = 2:N-1
        x = a + i * h;
        a_arr(i) = 2 + h^2 * q(x);
        b_arr(i) = -1 + (h/2) * p(x);
        c_arr(i) = -1 - (h/2) * p(x);
        d_arr(i) = -h^2 * r(x);
    end
    
    % Step 3: Last equation (i=N)
    x = b - h;
    a_arr(N) = 2 + h^2 * q(x);
    c_arr(N) = -1 - (h/2) * p(x);
    d_arr(N) = -h^2 * r(x) + (1 - (h/2) * p(x)) * beta;
    
    % Step 4: Solve tridiagonal system using Thomas algorithm
    l = zeros(N, 1);
    u = zeros(N, 1);
    z = zeros(N, 1);
    
    l(1) = a_arr(1);
    u(1) = b_arr(1) / l(1);
    z(1) = d_arr(1) / l(1);
    
    % Step 5: Forward elimination
    for i = 2:N-1
        l(i) = a_arr(i) - c_arr(i) * u(i-1);
        u(i) = b_arr(i) / l(i);
        z(i) = (d_arr(i) - c_arr(i) * z(i-1)) / l(i);
    end
    
    % Step 6: Last equation
    l(N) = a_arr(N) - c_arr(N) * u(N-1);
    z(N) = (d_arr(N) - c_arr(N) * z(N-1)) / l(N);
    
    % Step 7-8: Back substitution
    w = zeros(N+2, 1);
    w(1) = alpha;          % w0 = alpha
    w(N+2) = beta;         % w_{N+1} = beta
    w(N+1) = z(N);         % w_N = z_N
    
    for i = N-1:-1:1
        w(i+1) = z(i) - u(i) * w(i+2);
    end
    
    % Step 9: Create x-values vector
    x_vals = zeros(N+2, 1);
    for i = 0:N+1
        x_vals(i+1) = a + i * h;
    end
    
    % Return results
    x = x_vals;
end