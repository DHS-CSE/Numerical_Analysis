x = [0, 1, 3, 4, 6, 7, 9, 10];
fx = [8, 6, 5.5, 1.5, 1.1, 1, 0.9, 0.9];
dfx = [-2, -0.35, -0.67, -0.6, -0.04, -0.08, 0, 0];

% A polynomial of least degrree interpolating the value f(x) of given data
N = length(x);
A = zeros(N, N);

for i = 1:N
    for j = 1:N
        A(i, j) = x(i).^(j-1);
    end
end

coeff = A \ fx';

x_test = linspace(0, 10, 1001);
y_test = p_x(x_test);

figure
plot(x_test, y_test, 'b-', 'LineWidth', 2); hold on;
scatter(x, fx, 40, 'r', 'filled')
legend('Interpolating Polynomial', 'Data Points', 'Location', 'Best');
title('Interpolating Polynomial Fit to Given Data');
xlabel('x'); ylabel('Value');
grid on
set(gca, 'FontSize', 12);
hold off

%% A natural cubic spline that interpolates the value f(x) of given data
function sx = natural_cubic(x_data)
x = [0, 1, 3, 4, 6, 7, 9, 10];
fx = [8, 6, 5.5, 1.5, 1.1, 1, 0.9, 0.9];
N = length(x);
h = zeros(1, N-1);
for i = 1:N-1
    h(i) = x(i+1) - x(i);
end

A = zeros(N-2, N-2);
for i = 1:N-2
    A(i, i) = (h(i) + h(i+1)) / 3;
    if 1 < i
        A(i, i-1) = h(i-1) / 6;
    elseif i < N-2
        A(i, i+1) = h(i+1) / 6;
    end
end

D = zeros(N-2, 1);

for i = 1:N-2
    D(i) = (fx(i+2) - fx(i+1)) / h(i+1) - (fx(i+1) - fx(i)) / h(i);
end

M_inner = A \ D;
M = zeros(N, 1);
M(2:N-1) = M_inner;

N_data = length(x_data);    
sx = zeros(size(x_data));
for i = 1:N_data
    for j = 1:N-1
        if (x(j) <= x_data(i)) && (x_data(i) <= x(j+1))
            sx(i) = ((x(j+1) - x_data(i)).^3*M(j) + (x_data(i)-x(j)).^3*M(j+1))/(6*h(j)) + (fx(j)/h(j)-h(j)*M(j)/6)*(x(j+1)-x_data(i)) + (fx(j+1)/h(j)-h(j)*M(j+1)/6)*(x_data(i)-x(j));
        end
    end
end
end

cubic_approximation = natural_cubic(x_test);

figure
plot(x_test, cubic_approximation, 'b-', 'LineWidth', 2); hold on;
scatter(x, fx, 40, 'r', 'filled')
title('Natural Cubic Spline Interpolation')
xlabel('x'); ylabel('Value')
legend('Cubic Spline', 'Data Points', 'Location', 'Best')
grid on
set(gca, 'FontSize', 12)
hold off
%% A cubic Hermite that interpolates the value f(x) and the slope f'(x) of given data
function px = cubic_hermite(x_data)
x = [0, 1, 3, 4, 6, 7, 9, 10];
fx = [8, 6, 5.5, 1.5, 1.1, 1, 0.9, 0.9];
dfx = [-2, -0.35, -0.67, -0.6, -0.04, -0.08, 0, 0];
N = length(x);
N_data = length(x_data);

f_phi = zeros(size(x_data));
fprime_psi = zeros(size(x_data));
px = zeros(size(x_data));
for i = 1:N_data
    x_data_i = x_data(i);
    for j = 1:N-1
        if (x(j) <= x_data_i) && (x_data_i <= x(j+1))
            f_phi(i) = fx(j)*(x_data_i-x(j+1)).^2*((x(j)-x(j+1))+2*(x(j)-x_data_i))/(x(j)-x(j+1)).^3 + fx(j+1)*(x_data_i-x(j)).^2*((x(j+1)-x(j))+2*(x(j+1)-x_data_i))/(x(j+1)-x(j)).^3;
            fprime_psi(i) = dfx(j)*(x_data_i-x(j))*(x_data_i-x(j+1)).^2/(x(j)-x(j+1)).^2 + dfx(j+1)*(x_data_i-x(j)).^2*(x_data_i-x(j+1))/(x(j)-x(j+1)).^2;
        end
    end
    px(i) = f_phi(i) + fprime_psi(i);
end
end
cubic_hermite_approximation = cubic_hermite(x_test);

figure
plot(x_test, cubic_hermite_approximation, 'b-', 'LineWidth', 2); hold on;
scatter(x, fx, 40, 'r', 'filled');
title('Cubic Hermite Interpolation')
xlabel('x')
ylabel('Value')
legend('Cubic Hermite', 'Data Points', 'Location', 'Best')
grid on
set(gca, 'FontSize', 12)
hold off