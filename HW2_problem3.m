% Least squares approximation
for n = 1:4
    f = @(x) x.^(n+1);

    A = zeros(n+1, n+1);
    B = zeros(n+1, 1);

    for i = 0:n
        for j = 0:n
            A(i+1, j+1) = integral(@(x) x.^(i+j), 0, 3);
        end
        B(i+1) = integral(@(x) f(x).*x.^i, 0, 3);
    end

    c = A\B;
    
    p = @(x) polyval(flip(c'), x);

    x_vals = linspace(0, 3, 500);
    f_vals = f(x_vals);
    p_vals = p(x_vals);

    err_fun = @(x) (f(x) - p(x)).^2;
    L2_error = sqrt(integral(err_fun, 0, 3));

    figure(n);
    plot(x_vals, f_vals, 'b-', 'LineWidth', 2); hold on;
    plot(x_vals, p_vals, 'r--', 'LineWidth', 2);
    legend('f(x)', 'p^*(x)');
    title(['Least Squares Approximation for n = ', num2str(n)]);
    xlabel('x'); ylabel('Value');

    figure(n + 10);
    plot(x_vals, abs(f_vals - p_vals), 'k-', 'LineWidth', 2);
    title(['Error: |f(x) - p^*(x)|, n = ', num2str(n)]);
    xlabel('x'); ylabel('Error');

    fprintf('\n--- n = %d ---\n', n);
    fprintf('Least squares coefficients (low to high):\n');
    disp(c');
    fprintf('L2 Norm of Error: %.6f\n', L2_error);
end
%% Best uniform approximation using exchange algorithm
clear all;
close all;
clc;

a = 0;
b = 3;

n_values = [1, 2, 3, 4];

for n = n_values
    fprintf('\n======== n = %d ========\n', n);
    
    f = @(x) x.^(n+1);
    
    [p_coeff, error_norm, alternating_points] = remez_exchange_algorithm(f, n, a, b);
    
    % 결과 출력
    fprintf('Polynomial coefficients (High to Low): ');
    disp(p_coeff);
    fprintf('L∞ norm of error: %.10f\n', error_norm);
    fprintf('Alternating points: ');
    disp(alternating_points');

    plot_approximation(f, p_coeff, a, b, n);

    plot_error(f, p_coeff, a, b, n, alternating_points);
end

function [p_coeff, error_norm, alternating_points] = remez_exchange_algorithm(f, n, a, b, tol)

    if nargin < 5
        tol = 1e-12;
    end
    
    m = n + 2;
    cheb_nodes = cos(pi * ((0:m-1) + 0.5) / m)';
    x = ((b-a) * cheb_nodes + (b+a)) / 2;
    x = sort(x);

    max_iter = 100;
    iter = 0;
    prev_error = Inf;
    conv_count = 0;
    
    while iter < max_iter
        iter = iter + 1;

        [p_coeff, error] = compute_minimax_polynomial(f, x, n);

        p = @(t) polyval(p_coeff, t);

        error_func = @(t) f(t) - p(t);

        [max_errors, max_locs] = find_extrema(error_func, a, b, 100);
        
        [~, idx] = max(abs(max_errors));
        max_error = max_errors(idx);
        max_loc = max_locs(idx);

        if abs(abs(max_error) - abs(error)) < tol
            conv_count = conv_count + 1;
            if conv_count >= 3
                break;
            end
        else
            conv_count = 0;
        end

        x = update_exchange_points(x, max_loc, error_func);

        if abs(abs(error) - prev_error) < tol
            break;
        end
        prev_error = abs(error);
    end
    
    fprintf('Iteration: %d\n', iter);
    
    error_norm = abs(error);
    alternating_points = sort(x);
end

function [p_coeff, error] = compute_minimax_polynomial(f, x, n)
    
    m = length(x);
    A = zeros(m, n+2);
    b = zeros(m, 1);
    
    for i = 1:m
        for j = 1:n+1
            A(i, j) = x(i)^(n+1-j);
        end

        A(i, n+2) = (-1)^(i-1);
        b(i) = f(x(i));
    end
    
    z = A \ b;
    
    p_coeff = z(1:n+1)';
    error = z(n+2);
end

function [extrema_values, extrema_locs] = find_extrema(func, a, b, num_intervals)
    
    segment_length = (b - a) / num_intervals;
    intervals = a:segment_length:b;
    
    extrema_values = [];
    extrema_locs = [];
    
    for i = 1:length(intervals)-1
        start_point = intervals(i);
        end_point = intervals(i+1);
        
        % 극값 탐색
        options = optimset('Display', 'off', 'TolX', 1e-12);
        [x_max, fval] = fminbnd(@(x) -abs(func(x)), start_point, end_point, options);
        
        extrema_values = [extrema_values; func(x_max)];
        extrema_locs = [extrema_locs; x_max];
    end
    
    extrema_values = [extrema_values; func(a); func(b)];
    extrema_locs = [extrema_locs; a; b];
    
    [extrema_locs, unique_idx] = unique(extrema_locs);
    extrema_values = extrema_values(unique_idx);

    [~, sort_idx] = sort(abs(extrema_values), 'descend');
    extrema_values = extrema_values(sort_idx);
    extrema_locs = extrema_locs(sort_idx);
end

function new_x = update_exchange_points(x, new_point, error_func)
    
    candidate_x = sort([x; new_point]);
    
    n_plus_2 = length(x);
    
    best_x = select_best_alternating_set(candidate_x, error_func, n_plus_2);
    
    new_x = best_x;
end

function best_set = select_best_alternating_set(candidate_points, error_func, n_plus_2)
    
    errors = error_func(candidate_points);
    
    if length(candidate_points) <= n_plus_2
        best_set = candidate_points;
        return;
    end
    
    best_set = candidate_points;
    
    to_remove = length(best_set) - n_plus_2;
    
    for i = 1:to_remove
        best_set = remove_weakest_point(best_set, error_func);
    end
end

function reduced_set = remove_weakest_point(point_set, error_func)
    
    errors = error_func(point_set);
    m = length(point_set);
    
    same_sign_idx = [];
    for i = 2:m-1
        if sign(errors(i-1)) == sign(errors(i)) || sign(errors(i)) == sign(errors(i+1))
            same_sign_idx = [same_sign_idx; i];
        end
    end
  
    if length(same_sign_idx) == 0
        [~, min_idx] = min(abs(errors));
        reduced_set = point_set;
        reduced_set(min_idx) = [];
        return;
    end
    
    [~, min_idx] = min(abs(errors(same_sign_idx)));
    remove_idx = same_sign_idx(min_idx);
    
    reduced_set = point_set;
    reduced_set(remove_idx) = [];
end

function plot_approximation(f, p_coeff, a, b, n)
    figure;
    x = linspace(a, b, 1000);
    y_orig = f(x);
    y_approx = polyval(p_coeff, x);
    
    plot(x, y_orig, 'b-', 'LineWidth', 2);
    hold on;
    plot(x, y_approx, 'r--', 'LineWidth', 2);
    
    title(sprintf('Best Uniform Approximation (n = %d)', n));
    xlabel('x');
    ylabel('f(x)');
    legend('f(x) = x^{n+1}', 'Best Uniform Approximation');
    grid on;
end

function plot_error(f, p_coeff, a, b, n, alternating_points)
    figure;
    x = linspace(a, b, 1000);
    y_orig = f(x);
    y_approx = polyval(p_coeff, x);
    error = y_orig - y_approx;
    
    plot(x, abs(error), 'b-', 'LineWidth', 2);
    hold on;
    
    error_at_points = f(alternating_points) - polyval(p_coeff, alternating_points);
    plot(alternating_points, abs(error_at_points), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    title(sprintf('Error : |f(x) - p*(x)| (n = %d)', n));
    xlabel('x');
    ylabel('Error');
    grid on;
    
    max_error = max(abs(error));
    line([a, b], [max_error, max_error], 'Color', 'k', 'LineStyle', '--');
    
    legend('Error', 'Alternating points', 'Location', 'best');
    
    text(a + 1.5, max_error*0.8, sprintf('Max Error = %.6f', max_error), 'FontSize', 12, 'fontweight','bold', 'Color', 'r', 'HorizontalAlignment', 'center');
    
    fprintf('||f - p*||_\\infty = %.10f\n', max_error);
end
%% Least square approximation with f(x) = sin3x on [-pi, pi]
for n = 1:4
    f = @(x) sin(3*x);

    A = zeros(n+1, n+1);
    B = zeros(n+1, 1);

    for i = 0:n
        for j = 0:n
            A(i+1, j+1) = integral(@(x) x.^(i+j), -pi, pi);
        end
        B(i+1) = integral(@(x) f(x).*x.^i, -pi, pi);
    end

    c = A\B;
    
    p = @(x) polyval(flip(c'), x);

    x_vals = linspace(-pi, pi, 500);
    f_vals = f(x_vals);
    p_vals = p(x_vals);

    err_fun = @(x) (f(x) - p(x)).^2;
    L2_error = sqrt(integral(err_fun, 0, 3));

    figure(n);
    plot(x_vals, f_vals, 'b-', 'LineWidth', 2); hold on;
    plot(x_vals, p_vals, 'r--', 'LineWidth', 2);
    legend('f(x)', 'p^*(x)');
    title(['Least Squares Approximation for n = ', num2str(n)]);
    xlabel('x'); ylabel('Value');

    figure(n + 10);
    plot(x_vals, abs(f_vals - p_vals), 'k-', 'LineWidth', 2);
    title(['Error: |f(x) - p^*(x)|, n = ', num2str(n)]);
    xlabel('x'); ylabel('Error');

    fprintf('\n--- n = %d ---\n', n);
    fprintf('Least squares coefficients (low to high):\n');
    disp(c');
    fprintf('L2 Norm of Error: %.6f\n', L2_error);
end
%% Best uniform approximation using exchange algorithm(f(x) = sin3x)
clear all;
close all;
clc;

a = -pi;
b = pi;

n_values = [1, 2, 3, 4];

for n = n_values
    fprintf('\n======== n = %d ========\n', n);
    f = @(x) sin(3*x);

    [p_coeff, error_norm, alternating_points] = exchange_algorithm(f, n, a, b);
    
    fprintf('Polynomial coefficients (High to Low): ');
    disp(p_coeff);
    fprintf('L∞ norm of error: %.10f\n', error_norm);
    fprintf('Alternating Points: ');
    disp(alternating_points');

    plot_approximation_sin(f, p_coeff, a, b, n);

    plot_error_sin(f, p_coeff, a, b, n, alternating_points);
end

function [p_coeff, error_norm, alternating_points] = exchange_algorithm(f, n, a, b, tol)
    
    if nargin < 5
        tol = 1e-12;
    end
    
    m = n + 2;
    
    cheb_nodes = cos(pi * ((0:m-1) + 0.5) / m)';
    x = ((b-a) * cheb_nodes + (b+a)) / 2;
    x = sort(x);

    max_iter = 100;
    iter = 0;
    prev_error = Inf;
    conv_count = 0;
    
    while iter < max_iter
        iter = iter + 1;
        
        [p_coeff, error] = minimax_polynomial(f, x, n);
        
        p = @(t) polyval(p_coeff, t);

        error_func = @(t) f(t) - p(t);
        
        [max_errors, max_locs] = find_extrema_sin(error_func, a, b, 100);
        
        [~, idx] = max(abs(max_errors));
        max_error = max_errors(idx);
        max_loc = max_locs(idx);
        
        if abs(abs(max_error) - abs(error)) < tol
            conv_count = conv_count + 1;
            if conv_count >= 3
                break;
            end
        else
            conv_count = 0;
        end
        
        x = update_exchange_points_sin(x, max_loc, error_func);
        
        if abs(abs(error) - prev_error) < tol
            break;
        end
        prev_error = abs(error);
    end
    
    fprintf('Iteration: %d\n', iter);
    
    error_norm = abs(error);
    alternating_points = sort(x);
end

function [p_coeff, error] = minimax_polynomial(f, x, n)
    
    m = length(x);
    A = zeros(m, n+2);
    b = zeros(m, 1);
    
    for i = 1:m
        for j = 1:n+1
            A(i, j) = x(i)^(n+1-j);
        end

        A(i, n+2) = (-1)^(i-1);
        b(i) = f(x(i));
    end
    
    z = A \ b;
    
    p_coeff = z(1:n+1)';
    error = z(n+2);
end

function [extrema_values, extrema_locs] = find_extrema_sin(func, a, b, num_intervals)
    
    segment_length = (b - a) / num_intervals;
    intervals = a:segment_length:b;
    
    extrema_values = [];
    extrema_locs = [];
    
    for i = 1:length(intervals)-1
        start_point = intervals(i);
        end_point = intervals(i+1);
        
        options = optimset('Display', 'off', 'TolX', 1e-12);
        [x_max, fval] = fminbnd(@(x) -abs(func(x)), start_point, end_point, options);
        
        extrema_values = [extrema_values; func(x_max)];
        extrema_locs = [extrema_locs; x_max];
    end
    
    extrema_values = [extrema_values; func(a); func(b)];
    extrema_locs = [extrema_locs; a; b];
    
    [extrema_locs, unique_idx] = unique(extrema_locs);
    extrema_values = extrema_values(unique_idx);
    
    [~, sort_idx] = sort(abs(extrema_values), 'descend');
    extrema_values = extrema_values(sort_idx);
    extrema_locs = extrema_locs(sort_idx);
end

function new_x = update_exchange_points_sin(x, new_point, error_func)
    
    candidate_x = sort([x; new_point]);
    
    n_plus_2 = length(x);
    
    best_x = select_best_alternating_set_sin(candidate_x, error_func, n_plus_2);
    
    new_x = best_x;
end

function best_set = select_best_alternating_set_sin(candidate_points, error_func, n_plus_2)
    
    errors = error_func(candidate_points);
    
    if length(candidate_points) <= n_plus_2
        best_set = candidate_points;
        return;
    end
    
    best_set = candidate_points;
    
    to_remove = length(best_set) - n_plus_2;
    
    for i = 1:to_remove
        best_set = remove_weakest_point_sin(best_set, error_func);
    end
end

function reduced_set = remove_weakest_point_sin(point_set, error_func)
    
    errors = error_func(point_set);
    m = length(point_set);
    
    same_sign_idx = [];
    for i = 2:m-1
        if sign(errors(i-1)) == sign(errors(i)) || sign(errors(i)) == sign(errors(i+1))
            same_sign_idx = [same_sign_idx; i];
        end
    end
    
    if length(same_sign_idx) == 0

        [~, min_idx] = min(abs(errors));
        reduced_set = point_set;
        reduced_set(min_idx) = [];
        return;
    end
    
    [~, min_idx] = min(abs(errors(same_sign_idx)));
    remove_idx = same_sign_idx(min_idx);
    
    reduced_set = point_set;
    reduced_set(remove_idx) = [];
end

function plot_approximation_sin(f, p_coeff, a, b, n)
    figure;
    x = linspace(a, b, 1000);
    y_orig = f(x);
    y_approx = polyval(p_coeff, x);
    
    plot(x, y_orig, 'b-', 'LineWidth', 2);
    hold on;
    plot(x, y_approx, 'r--', 'LineWidth', 2);
    
    title(sprintf('Best Uniform Approximation: (n = %d)', n));
    xlabel('x');
    ylabel('f(x)');
    legend('f(x) = sin3x', 'Best Uniform Approximation');
    grid on;
end

function plot_error_sin(f, p_coeff, a, b, n, alternating_points)

    figure;
    x = linspace(a, b, 1000);
    y_orig = f(x);
    y_approx = polyval(p_coeff, x);
    error = y_orig - y_approx;
    
    plot(x, abs(error), 'b-', 'LineWidth', 2);
    hold on;
    
    error_at_points = f(alternating_points) - polyval(p_coeff, alternating_points);
    plot(alternating_points, abs(error_at_points), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    title(sprintf('Error :  |f(x) - p*(x)| (n = %d)', n));
    xlabel('x');
    ylabel('Error');
    grid on;
    
    max_error = max(abs(error));
    line([a, b], [max_error, max_error], 'Color', 'k', 'LineStyle', '--');
    
    legend('Error', 'Alternating Points', 'Location', 'best');
    
    text(a + 3, max_error*0.8, sprintf('Max Error = %.6f', max_error), 'FontSize', 12, 'fontweight','bold', 'Color', 'r', 'HorizontalAlignment', 'center');
    
    fprintf('||f - p*||_\\infty = %.10f\n', max_error);
end