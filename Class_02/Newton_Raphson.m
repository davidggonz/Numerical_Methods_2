%Newton-Raphson Method
%f'(x0) = (f(xn) - f(xn+1))/(xn - xn+1)

%If f(x) with x = xn intercepts the x-axis in xn+1
%f'(xn) = f(xn)/(xn - xn+1)

%For n equations 
%xn+1 = xn - f(x)/f'(x)
% --> xn+1 = xn + v, in which Jv = -f(xn)
%J being the Jacobian of the system of equations
%So v = -J^{-1}f(xn)
%Which leads to: xn+1 = 

%Examples:

% F = ( 4x^2 + y^2 + 2xy - y - 2 --> J ( 8x+2y  2y+2x-1 )
%       2x^2 + 3xy + y^2 - 3 )         (  4x+3y  3x+2y   )

x = 0.1;
y = 0.1;

solution_n = [x,y]';
tol = 1e-10;
N = 100;

for i=1:N
    J = [8*x + 2*y, 2*y+2*x-1;
        4*x+3*y,3*x+2*y];
    F = [4*x^2+y^2+2*x*y-y-2;
        2*x^2+3*x*y+y^2-3];
    v = J\(-F);
    snn = solution_n + v;
    test = sqrt(dot(solution_n - snn, solution_n - snn));
    if test <= tol
        break;
    end
    solution_n = snn;
    x = solution_n(1);
    y = solution_n(2);
end

%Being U a continous function with continous derivatives
%By Taylor:
%U(x+h) = U(x) +hU'(x) +1/2 h^2U''(x) +1/6h^3U'''(x) +...
%U(x-h) = U(x) -hU'(x) +1/2 h^2U''(x) -1/6h^3U'''(x) +...
%So U(x+h)+U(x-h) = 2U(x) + h^2U''(x) + O(h^4)

%U''(x) = 1/h^2 [U(x+h) -2U(x) + U(x-h)], with and error of h^
%So centered: U'(x) = 1/2h [U(x+h) - U(x-h)], with and error of h

%Being y" = P(x)y' + q(x)y + r(x)
%For a =< x =< b y(a) = α & y(b) = β
