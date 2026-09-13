%% Derve Gauss-Kronrod pair, Verify what the highest degree of polynommial is
syms x a b c d e
P = a*x^4 + b*x^3 + c*x^2 + d*x + e;
W = 4*x^3 - 3*x;
weight = 1/sqrt(1 - x^2);

eqs = sym(zeros(5,1));

for k = 0:3
    integrand = x^k * W * P * weight;
    eqs(k+1) = int(integrand, x, -1, 1);
end

vars = [a, b, c, d, e];
[A, B] = equationsToMatrix(eqs, vars);

nullspace = null(double(A));

kronrod_points = roots(nullspace.');
kronrod_points = kronrod_points.';

n = 3;
k = 1:n;
nodes_gauss = cos((2*k - 1) * pi / (2*n));
nodes = [nodes_gauss kronrod_points];
disp('Gauss_points_3:')
disp(nodes_gauss)
disp('Kronrod_points_4:')
disp(kronrod_points)
disp('Entire_points:')
disp(nodes)

syms x
w = 1 ./ (1-x.^2).^(1/2);

W = zeros(2*n + 1, 1);

for i=1:(2*n + 1)
    bottom_term = 1;
    for j = 1:(2*n + 1)
        if i~=j
            bottom_term = bottom_term * (nodes(i) - nodes(j));
        end
    end
    coefficients = poly(nodes([1:i-1, i+1:(2*n + 1)])) / bottom_term;
    weight_function = w * poly2sym(coefficients, x);
    W(i) = int(weight_function, x , -1, 1);
end

tol = 1e-8;
k = 0;
while true
    syms x
    y = x.^k ./ (1-x.^2).^(1/2);
    exact = int(y, x, -1, 1);

    approximation = W.' * nodes.^k.';
    err = abs(approximation - exact);
    fprintf('%d\t %.4e\n', k, err);

    if err < tol
        k = k + 1;
    else
        max_degree = k-1;
        break;
    end
end

fprintf('\n The highest degree of polynomial: %d\n', max_degree);
%% Find the qpproximation of integral by using Gauss-Kronrod
syms x a b c d e
P = a*x^4 + b*x^3 + c*x^2 + d*x + e;
W = 5*x^3 - 3*x;

eqs = sym(zeros(5,1));

for k = 0:3
    integrand = x^k * W * P;
    eqs(k+1) = int(integrand, x, -1, 1);
end

vars = [a, b, c, d, e];
[A, B] = equationsToMatrix(eqs, vars);

nullspace = null(double(A));

kronrod_points = roots(nullspace.');
kronrod_points = kronrod_points.';
n_kr = length(kronrod_points);

n = 3;
nodes_gauss = roots([5, 0, -3, 0]).';
nodes = [nodes_gauss kronrod_points];

W = zeros(2*n + 1, 1);

for i=1:(2*n + 1)
    bottom_term = 1;
    for j = 1:(2*n + 1)
        if i~=j
            bottom_term = bottom_term * (nodes(i) - nodes(j));
        end
    end
    coefficients = poly(nodes([1:i-1, i+1:(2*n + 1)])) / bottom_term;
    weight_function = poly2sym(coefficients, x);
    W(i) = int(weight_function, x , -1, 1);
end

a = 0;
b = pi/4;
f = @(t) exp( (b-a)/2 * t + (b+a)/2 ) .* cos( 4 * ( (b-a)/2 * t + (b+a)/2 ) ) ./ sqrt(1 - ( (b-a)/2 * t + (b+a)/2 ).^2 );
I = (b-a) / 2 * W.' * f(nodes.');
fprintf('The approximation of integral: %.8f\n', double(I));

syms x 
f_value =  exp(x) .* ...
     cos(4*x) ./ ...
     sqrt(1 - ( x ).^2 );
exact_value = int(f_value, x, 0, pi/4);
fprintf('Exact: %.8f\n', double(exact_value));

fprintf('Error: %d\n', abs(I - double(exact_value)))