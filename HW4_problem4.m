function shooting_method_secant
    % Step size
    h = 1/16;
    t = 0:h:1;
    N = length(t);

    % Tolerance
    tol = 1e-3;

    % Initial guesses for y(0)
    s0 = -1;
    s1 = 1;

    % Secant iteration
    max_iter = 50;
    for iter = 1:max_iter
        F0 = boundary_mismatch(s0, t, h);
        F1 = boundary_mismatch(s1, t, h);

        % Secant update
        s_new = s1 - F1*(s1 - s0)/(F1 - F0);

        fprintf("Iteration %d: s = %.6f, boundary error = %.6f\n", iter, s_new, F1);

        % Check for convergence
        if abs(s_new - s1) < tol
            break;
        end

        s0 = s1;
        s1 = s_new;
    end

    % Final integration with best initial value
    [~, y1_final] = integrate_ode(s_new, t, h);

    % Plot result
    plot(t, y1_final, '-o');
    xlabel('t'); ylabel('y(t)');
    title('Solution using Shooting Method with Secant');
    grid on;
end

function mismatch = boundary_mismatch(s, t, h)
    [y2, y1] = integrate_ode(s, t, h);
    % y(1) + y'(1) should be 3
    mismatch = y1(end) + y2(end) - 3;
end

function [y2, y1] = integrate_ode(s, t, h)
    N = length(t);
    y1 = zeros(1, N);  % y
    y2 = zeros(1, N);  % y'

    % Initial conditions
    y1(1) = s;
    y2(1) = 0;

    for i = 1:N-1
        % Define f1 and f2
        f1 = y2(i);
        if y1(i) == 0
            f2 = 0; % Avoid division by zero
        else
            f2 = ((y2(i))^2 + 2*(t(i))^2 * y1(i)) / (2*y1(i));
        end

        % Euler method
        y1(i+1) = y1(i) + h * f1;
        y2(i+1) = y2(i) + h * f2;
    end
end

shooting_method_secant
%% finite difference method
function finite_difference_nonlinear_bvp
    ns = [8, 16, 32]; % Number of subintervals
    
    figure;
    hold on;
    colors = ['r', 'g', 'b'];
    
    for idx = 1:length(ns)
        n = ns(idx);
        h = 1 / n;
        t = linspace(0, 1, n+1)';
        
        % Initial guess: linear guess satisfying roughly the BCs
        y_init = linspace(1, 2, n+1)';
        
        % Set options for fsolve correctly
        options = optimoptions('fsolve', ...
                               'Display', 'off', ...
                               'FunctionTolerance', 1e-8, ...
                               'StepTolerance', 1e-8);
                           
        % Solve nonlinear system using fsolve
        y_sol = fsolve(@(y) nonlinear_system(y, t, h, n), y_init, options);
        
        plot(t, y_sol, [colors(idx) '-o'], 'DisplayName', sprintf('n = %d', n));
    end
    
    xlabel('t');
    ylabel('y(t)');
    title('Finite Difference Solution of Nonlinear BVP');
    legend show;
    grid on;
    hold off;
end

function F = nonlinear_system(y, t, h, n)
    % Nonlinear system based on finite difference discretization
    
    F = zeros(n+1, 1);
    
    % Boundary condition at t=0: y'(0) = 0 approximated by forward difference
    F(1) = (y(2) - y(1)) / h; % should be zero
    
    % Interior points i=2,...,n
    for i = 2:n
        y_i = y(i);
        y_ip1 = y(i+1);
        y_im1 = y(i-1);
        t_i = t(i);
        
        y_prime = (y_ip1 - y_im1) / (2*h);
        y_double_prime = (y_ip1 - 2*y_i + y_im1) / (h^2);
        
        F(i) = - y_prime^2 - 2*(t_i)^2*y_i + 2*y_i*y_double_prime;
    end
    
    % Boundary condition at t=1: y(1) + y'(1) = 3 approximated by backward difference
    y_prime_end = (y(n+1) - y(n)) / h;
    F(n+1) = y(n+1) + y_prime_end - 3;
end

finite_difference_nonlinear_bvp
%% collocation point method with cubic B-spline
function collocation_bvp_bspline
    clf;
    ns = [8, 16, 32];
    colors = ['r', 'g', 'b'];
    hold on;

    for idx = 1:length(ns)
        n = ns(idx);
        [t_vals, y_vals] = collocation_method(n);
        plot(t_vals, y_vals, 'Color', colors(idx), 'DisplayName', ['n = ', num2str(n)]);
    end

    title('Collocation Method Using Cubic B-splines');
    xlabel('t');
    ylabel('y(t)');
    legend;
    grid on;
    hold off;
end

function [t_vals, y_vals] = collocation_method(n)
    k = 4;  % Cubic B-spline: degree = k - 1
    num_basis = n + k - 2;
    num_knots = num_basis + k;

    % Uniform knot vector (clamped)
    knots = [zeros(1, k), linspace(0, 1, num_knots - 2*k), ones(1, k)];

    % Greville abscissae as collocation points
    greville = zeros(1, num_basis);
    for i = 1:num_basis
        greville(i) = sum(knots(i+1:i+k-1)) / (k - 1);
    end

    % Initial guess for coefficients
    c0 = ones(num_basis, 1);

    % Solve nonlinear system via fsolve
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    c_sol = fsolve(@(c) residuals(c, knots, greville, k), c0, options);

    % Evaluate final spline
    t_vals = linspace(0, 1, 200);
    y_vals = zeros(size(t_vals));
    for i = 1:num_basis
        y_vals = y_vals + c_sol(i) * bspline_basis(i, k, knots, t_vals);
    end
end

function R = residuals(c, knots, t_colloc, k)
    num_basis = length(c);
    y = zeros(size(t_colloc));
    yp = zeros(size(t_colloc));
    ypp = zeros(size(t_colloc));

    for i = 1:num_basis
        B = bspline_basis(i, k, knots, t_colloc);
        dB = bspline_derivative(i, k, knots, t_colloc, 1);
        ddB = bspline_derivative(i, k, knots, t_colloc, 2);

        y = y + c(i) * B;
        yp = yp + c(i) * dB;
        ypp = ypp + c(i) * ddB;
    end

    % Residuals from the differential equation
    R = -(yp.^2) - 2.*(t_colloc.^2).*y + 2.*y.*ypp;

    % Enforce boundary conditions
    % y'(0) = 0
    d0 = 0;
    for i = 1:num_basis
        d0 = d0 + c(i) * bspline_derivative(i, k, knots, 0, 1);
    end
    R(end+1) = d0;

    % y(1) + y'(1) = 3
    y1 = 0;
    dy1 = 0;
    for i = 1:num_basis
        y1 = y1 + c(i) * bspline_basis(i, k, knots, 1);
        dy1 = dy1 + c(i) * bspline_derivative(i, k, knots, 1, 1);
    end
    R(end+1) = y1 + dy1 - 3;
end

function B = bspline_basis(i, k, knots, x)
    B = bspline_basis_rec(i, k, knots, x);
end

function dB = bspline_derivative(i, k, knots, x, order)
    B = bspline_basis_rec(i, k, knots, x);
    for d = 1:order
        B = diff_bspline(i, k - d + 1, knots, x);
    end
    dB = B;
end

function B = bspline_basis_rec(i, k, t, x)
    if k == 1
        B = double(t(i) <= x & x < t(i+1));
        if t(i+1) == 1  % right end inclusion
            B(x == 1) = 1;
        end
    else
        denom1 = t(i+k-1) - t(i);
        denom2 = t(i+k) - t(i+1);

        term1 = 0;
        term2 = 0;
        if denom1 > 0
            term1 = ((x - t(i)) ./ denom1) .* bspline_basis_rec(i, k-1, t, x);
        end
        if denom2 > 0
            term2 = ((t(i+k) - x) ./ denom2) .* bspline_basis_rec(i+1, k-1, t, x);
        end
        B = term1 + term2;
    end
end

function dB = diff_bspline(i, k, t, x)
    denom1 = t(i+k-1) - t(i);
    denom2 = t(i+k) - t(i+1);

    term1 = 0;
    term2 = 0;
    if denom1 > 0
        term1 = (k-1)/denom1 * bspline_basis_rec(i, k-1, t, x);
    end
    if denom2 > 0
        term2 = -(k-1)/denom2 * bspline_basis_rec(i+1, k-1, t, x);
    end
    dB = term1 + term2;
end

collocation_bvp_bspline