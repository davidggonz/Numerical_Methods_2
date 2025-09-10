function [w, iter] = newton_raphson_nonlinear(F_func, J_func, w0, h, alpha, beta, tol, max_iter)
    % Newton-Raphson method for nonlinear finite difference equations
    w = w0(:); 
    N = length(w);
    iter = 0;
    
    for k = 1:max_iter
        iter = k;
        F = F_func(w, h, alpha, beta);
        J = J_func(w, h, alpha, beta);
        
        % Solve J * v = -F
        v = J \ (-F);
        
        w_new = w + v;
        
        % Check convergence
        if norm(w_new - w, inf) < tol
            w = w_new;
            fprintf('Converged in %d iterations\n', iter);
            return;
        end
        
        w = w_new;
    end
    
    fprintf('Did not converge after %d iterations\n', max_iter);
end