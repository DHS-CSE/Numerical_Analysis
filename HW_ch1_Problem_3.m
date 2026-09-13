% Problem 3
% Define each g_i(x) function
g1 = @(x) (x^2 + 2) / 3;
g2 = @(x) sqrt((3 * x - 2));
g3 = @(x) 3 - 2 / x;
g4 = @(x) (x^2 - 2) / (2 * x - 3);

% Initial guess and parameters
x0 = 2.5;

% Fixed-point iteration function
function [x, n, convergence_rate] = fixed_point(g, x0, tol, max_iter)
    x = x0;
    convergence_rate = [];
    for n = 1:max_iter
        x_new = g(x);

        convergence_value = log10(x_new - 2) / log10(x - 2);
        convergence_rate = [convergence_rate; convergence_value];

        if abs(x_new - x) < tol
            x = x_new;
            return;
        end
        x = x_new;
    end
end

% Run the fixed-point iteration for each g_i
[x1, n1, convergence_rate1] = fixed_point(g1, x0, tolerance, max_iter);
fprintf('g1: x = %.6f, iterations = %d\n', x1, n1);

[x2, n2, convergence_rate2] = fixed_point(g2, x0, tolerance, max_iter);
fprintf('g2: x = %.6f, iterations = %d\n', x2, n2);
disp('convergence_rate of g2 : ')
disp(convergence_rate2)

[x3, n3, convergence_rate3] = fixed_point(g3, x0, tolerance, max_iter);
fprintf('g3: x = %.6f, iterations = %d\n', x3, n3);
disp('convergence_rate of g3 : ')
disp(convergence_rate3)

[x4, n4, convergence_rate4] = fixed_point(g4, x0, tolerance, max_iter);
fprintf('g4: x = %.6f, iterations = %d\n', x4, n4);
disp('convergence_rate of g4 : ')
disp(convergence_rate4)