%% adaptive Simpson's quadrature rule
syms x
y = 4 / (1 + x^2);
exact_int = double(int(y, x, 0, 2));
fprintf('Exact_value = %.10f\n', exact_int);
f = matlabFunction(y);
tol = 1e-3;

a = 0;
b = 2;

function s = simpson(f, a, b)
    c = (a + b)/2;
    s = (b - a)/6 * (f(a) + 4*f(c) + f(b));
end

function [result, intervals, parts] = adaptive_simpson(f, a, b, tol, S, depth)
    c = (a + b)/2;
    S_left = simpson(f, a, c);
    S_right = simpson(f, c, b);

    if abs((S_left + S_right)-S)/15 < tol

        result = S_left + S_right + (S_left + S_right - S)/15;
        intervals = 2;
        parts = [S_left, S_right];
    else
        [res1, int1, parts1] = adaptive_simpson(f, a, c, tol/2, S_left, depth+1);
        [res2, int2, parts2] = adaptive_simpson(f, c, b, tol/2, S_right, depth+1);
        result = res1 + res2;
        intervals = int1 + int2;
        parts = [parts1, parts2];
    end
end

initial_S = simpson(f, a, b);
[approximation, final_intervals, partial_sums] = adaptive_simpson(f, a, b, tol, initial_S, 0);
error = abs(approximation - exact_int);

fprintf('Approximation = %.10f\n', approximation);
fprintf('Error = %.10f\n', error);
fprintf('Total final intervals = %d\n', final_intervals);
fprintf('Partial sums per interval:\n');

disp(partial_sums')
%% midpoint rule and trapezoid rule
clear;
close all;
clc;

syms x
y = 4 ./ (1 + x.^2);
exact_int = int(y, x, 0, 2);
exact_int = double(exact_int);
fprintf('Exact_value = %.10f\n', exact_int);
f = matlabFunction(y);

tol = 1e-3;

function [Approximation_mid, Approximation_tr, Partial_sum_mid, Partial_sum_tr, Intervals] = adaptive_mid_tr(f, a, b, tol)
    c = (a+b)/2;
    Partial_sum_mid = (b-a)*f(c);
    Partial_sum_tr = (b-a)/2*(f(a)+f(b));

    if abs(Partial_sum_tr - Partial_sum_mid) > tol
        [Approximation_mid1, Approximation_tr1, Partial_sum_mid1, Partial_sum_tr1, Int1] = adaptive_mid_tr(f, a, c, tol/2);
        [Approximation_mid2, Approximation_tr2, Partial_sum_mid2, Partial_sum_tr2, Int2] = adaptive_mid_tr(f, c, b, tol/2);
        Intervals = Int1 + Int2;
        Approximation_mid = Approximation_mid1 + Approximation_mid2;
        Approximation_tr = Approximation_tr1 + Approximation_tr2;
        Partial_sum_mid = [Partial_sum_mid1 Partial_sum_mid2];
        Partial_sum_tr = [Partial_sum_tr1 Partial_sum_tr2];
    else
        Intervals = 1;
        Approximation_mid = Partial_sum_mid;
        Approximation_tr = Partial_sum_tr;
    end
end

[Approximation_mid, Approximation_tr, Partial_sum_mid, Partial_sum_tr, Intervals] = adaptive_mid_tr(f, 0, 2, tol);

err_mid = abs(exact_int - Approximation_mid);
err_tr = abs(exact_int - Approximation_tr);
fprintf('Interval = %d\n', Intervals)
fprintf('Error_mid = %d\n', err_mid)
fprintf('Error_tr = %d\n', err_tr)
fprintf('Approximation_mid = %d\n', Approximation_mid)
fprintf('Approximation_tr = %d\n', Approximation_tr)
fprintf('Partial sums on midpoint rule per interval:\n');
disp(Partial_sum_mid.')
fprintf('Partial sums on trapezoid rule per interval:\n');
disp(Partial_sum_tr.')