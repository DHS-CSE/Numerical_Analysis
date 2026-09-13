% Define functions
f1 = @(x) (x(1) + 3) * (x(2)^3 - 7) + 18;
f2 = @(x) sin(x(2) * exp(x(1)) - 1);
F = @(x) [f1(x); f2(x)];

% Jacobian matrix
df1dx1 = @(x) x(2)^3 - 7;
df1dx2 = @(x) 3 * (x(1) + 3) * x(2)^2;
df2dx1 = @(x) cos(x(2)*exp(x(1)) - 1) * x(2) * exp(x(1));
df2dx2 = @(x) cos(x(2)*exp(x(1)) - 1) * exp(x(1));
Jacobian = @(x) [df1dx1(x), df1dx2(x);
           df2dx1(x), df2dx2(x)];

% Newton method
x0 = [-0.2 1.4];
tolerance = 1e-6;
max_iter = 1000;
x_newton = x0;
convergence_rate_newton = zeros(max_iter, 1);
x_exact = [0 1];

for k = 1:max_iter
     F_matrix = F(x_newton);

     J_matrix = Jacobian(x_newton);

     delta_x = J_matrix \ (-F_matrix);

     x_k1_newton = x_newton + delta_x.';

     if norm(x_k1_newton - x_newton) < tolerance
        difference_k1 = x_k1_newton - x_exact;
        difference_k = x_newton - x_exact;
        convergence_value = log10(norm(difference_k1)) / log10(norm(difference_k));
        convergence_rate_newton(k, :) = convergence_value;
        convergence_rate_newton = convergence_rate_newton(1:k, :);
        x_newton_pred = x_k1_newton;
        Newton_iter = k;
        break;
     end

     difference_k1 = x_k1_newton - x_exact;
     difference_k = x_newton - x_exact;
     convergence_value = log10(norm(difference_k1)) / log10(norm(difference_k));
     convergence_rate_newton(k, :) = convergence_value;
     x_newton = x_k1_newton;

end

disp('Newton_method :')
disp(x_newton_pred);
disp('convergence_rate : ')
disp(convergence_rate_newton)
disp('Iteration of Newton Method : ')
disp(Newton_iter)

% Broyden method
x_broyden = x0;

H0 = Jacobian(x_broyden);
B_k = H0;
convergence_rate_broyden = [];

for k = 1:max_iter
    F_matrix = F(x_broyden);
    delta_x = B_k \ (-F_matrix);
    x_k1_broyden = x_broyden + delta_x.';

    difference_k1 = x_k1_broyden - x_exact;
    difference_k = x_broyden - x_exact;
    convergence_value = log10(norm(difference_k1)) / log10(norm(difference_k));
    convergence_rate_broyden = [convergence_rate_broyden; convergence_value];

    if norm(x_k1_broyden - x_broyden, inf) < tolerance
        x_broyden_pred =x_k1_broyden;
        convergence_rate_broyden = convergence_rate_broyden(1:k, :);
        Brodyen_iter = k;
        break;
    end

    s_k = x_k1_broyden - x_broyden;
    f1_k1 = (x_k1_broyden(1,1) +3) * (x_k1_broyden(1,2)^3 - 7) + 18;
    f1_k = (x_broyden(1,1) +3) * (x_broyden(1,2)^3 - 7) + 18;
    f2_k1 = sin(x_k1_broyden(1,2)*exp(x_k1_broyden(1,1)) - 1);
    f2_k = sin(x_broyden(1,2)*exp(x_broyden(1,1)) - 1);
    y_k = [(f1_k1 - f1_k);
        (f2_k1 - f2_k)];

    B_k1 = B_k + ((y_k - B_k * s_k.')*s_k) / (s_k * s_k.');

    B_k = B_k1;

    x_broyden = x_k1_broyden;
end
    
disp('Broyden_method :')
disp(x_broyden_pred);
disp('convergence_rate : ')
disp(convergence_rate_broyden)
disp('Iteration of Broyden method')
disp(Brodyen_iter)

% SOR-Newton method

% x = x0;
% omega = 1.5;
% 
% for iter = 1:max_iter
%     x_old = x;
%     x_gs = x;
% 
%     % === Inner loop: Gauss-Seidel Newton ===
%     for inner = 1:max_iter
%         x_prev = x_gs;
% 
%         % --- Update x1 using x2 ---
%         f1 = @(x1) (x1 + 3)*(x_gs(2)^2 - 7) + 18;
%         df1 = @(x1) (x_gs(2)^2 - 7);
%         x1_new = x_gs(1) - f1(x_gs(1)) / df1(x_gs(1));
% 
%         % --- Update x2 using updated x1 ---
%         f2 = @(x2) sin(x2 * exp(x1_new) - 1);
%         df2 = @(x2) cos(x2 * exp(x1_new) - 1) * exp(x1_new);
%         x2_new = x_gs(2) - f2(x_gs(2)) / df2(x_gs(2));
% 
%         x_gs = [x1_new; x2_new];
% 
%         % Check convergence of inner Newton loop
%         if norm(x_gs - x_prev) < tolerance
%             break;
%         end
%     end
% 
%     % === SOR Update ===
%     x = x_old + omega * (x_gs - x_old);
% 
%     if norm(x- x_old) < tolerance
%         fprintf('\nConverged in %d iterations.\n', iter);
%         fprintf('Final x = [%.6f %.6f]\n', x(1), x(2));
%         return;
%     end
% end
% 
% disp('Maximum iterations reached without convergence.');


x = x0;
omega = 1.5;

for iter = 1:max_iter
    x_old = x;
    x_gs = x;

    % === Inner loop: Gauss-Seidel Newton ===
    x_prev = x_gs;

    % --- Update x1 using x2 ---
    f1 = @(x1) (x1 + 3)*(x_gs(2)^2 - 7) + 18;
    df1 = @(x1) (x_gs(2)^2 - 7);
    x1_new = x_gs(1) - f1(x_gs(1)) / df1(x_gs(1));

    % --- Update x2 using updated x1 ---
    f2 = @(x2) sin(x2 * exp(x1_new) - 1);
    df2 = @(x2) cos(x2 * exp(x1_new) - 1) * exp(x1_new);
    x2_new = x_gs(2) - f2(x_gs(2)) / df2(x_gs(2));

    x_gs = [x1_new; x2_new];

    % === SOR Update ===
    x = x_old + omega * (x_gs - x_old);

    if norm(x- x_old) < tolerance
        fprintf('\nConverged in %d iterations.\n', iter);
        fprintf('Final x = [%.6f %.6f]\n', x(1), x(2));
        return;
    end
end