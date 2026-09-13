f =  @(x) 1 ./ (1 + x.^2);
x_10 = linspace(-5, 5, 11);
y_10 = f(x_10);
N_10 = length(y_10);
x_16 = linspace(-5, 5, 17);
y_16 = f(x_16);
N_16 = length(y_16);
x_test = linspace(-5, 5, 1001);
y_test = f(x_test);

% (a) A lagrange interpolation polynomial of degree 10 and degreee 16

L_10 = zeros(N_10, N_10);
for i=1:N_10
    bottom_term = 1;
    for j = 1:N_10
        if i~=j
            bottom_term = bottom_term * (x_10(i) - x_10(j));
        end
    end
    V_10 = poly(x_10([1:i-1, i+1:N_10])) / bottom_term;
    L_10(i, :) = y_10(i) * V_10;
end
P_10 = sum(L_10, 1);

y_approximation_10 = polyval(P_10, x_test);

L_16 = zeros(N_16, N_16);
for i=1:N_16
    bottom_term = 1;
    for j = 1:N_16
        if i~=j
            bottom_term = bottom_term * (x_16(i) - x_16(j));
        end
    end
    V_16 = poly(x_16([1:i-1, i+1:N_16])) / bottom_term;
    L_16(i, :) = y_16(i) * V_16;
end
P_16 = sum(L_16, 1);

y_approximation_16 = polyval(P_16, x_test);

figure(1)
plot(x_test, y_test, 'k-', 'Linewidth', 2); hold on;
plot(x_test, y_approximation_10, 'r--', 'Linewidth', 2); hold on;
legend('f(x)', 'p_{10}(x)');
xlabel('x'); ylabel('Value');
title('Runge Function and Lagrange Interpolation (Degree 10)');
grid on; hold off

figure(2)
plot(x_test, y_test, 'k-', 'Linewidth', 2); hold on;
plot(x_test, y_approximation_16, 'b--', 'Linewidth', 2);
legend('f(x)', 'p_{16}(x)');
xlabel('x'); ylabel('Value');
title('Runge Function and Lagrange Interpolation (Degree 16)');
grid on; hold off

error_10 = y_test - y_approximation_10;
error_16 = y_test - y_approximation_16;

max_error_10 = max(abs(error_10));
max_error_16 = max(abs(error_16));
x_center = mean(x_test);

figure(11)
plot(x_test, error_10, 'k-', 'Linewidth', 2); hold on;
title('Error: f(x) - p_{10}(x)');
xlabel('x'); ylabel('Error');
text(x_center, max(error_10)*0.9, sprintf('Max Error = %.4e', max_error_10), 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
grid on; hold off

figure(12)
plot(x_test, error_16, 'k-', 'Linewidth', 2); hold on;
title('Error: f(x) - p_{16}(x)');
xlabel('x'); ylabel('Error');
text(x_center, max(error_16)*0.9, sprintf('Max Error = %.4e', max_error_16), ...
    'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
grid on; hold off
%% (b) An interpolating polynomial of degree 10 and of degree 16
% using Chebyshev points

N_10 = 10;
x_cheby_10 = cos((2*(0:N_10) + 1)*pi / (2*N_10 + 2));
f_transfer =  @(x) 1 ./ (1 + (5*x).^2);
y_cheby_10 = f_transfer(x_cheby_10);

A = zeros(N_10+1, N_10+1);
for i = 0:N_10
    for j = 0:N_10
        A(i+1,j+1) = x_cheby_10(i+1).^j;
    end
end

coeff = A \ y_cheby_10';

p_cheby_10 = @(x) polyval(flip(coeff), x);

x_test = linspace(-1, 1, 1001);
p_10 = p_cheby_10(x_test);
f_test = f_transfer(x_test);

figure(1)
plot(x_test, p_10)
hold off
figure(2)
plot(x_test, f_test)
hold off

N_16 = 16;
x_cheby_16 = cos((2*(0:N_16) + 1)*pi / (2*N_16 + 2));
y_cheby_16 = f_transfer(x_cheby_16);

A = zeros(N_16+1, N_16+1);
for i = 0:N_16
    for j = 0:N_16
        A(i+1,j+1) = x_cheby_16(i+1).^j;
    end
end

coeff = A \ y_cheby_16';

p_cheby_16 = @(x) polyval(flip(coeff), x);

p_16 = p_cheby_16(x_test);

figure(1)
plot(x_test, f_test, 'k-', 'LineWidth', 2); hold on;
plot(x_test, p_10, 'r--', 'LineWidth', 2);
legend('f(x)', 'p_{10}(x)');
title('f(x) and p_{10}(x) with Chebyshev Nodes');
xlabel('x'); ylabel('Value');
grid on; hold off

figure(2)
plot(x_test, f_test, 'k-', 'LineWidth', 2); hold on;
plot(x_test, p_16, 'b--', 'LineWidth', 2);
legend('f(x)', 'p_{16}(x)');
title('f(x) and p_{16}(x) with Chebyshev Nodes');
xlabel('x'); ylabel('Value');
grid on; hold off

% Error plots
error_10 = f_test - p_10;
error_16 = f_test - p_16;

max_error_10 = max(abs(error_10));
max_error_16 = max(abs(error_16));
x_center = mean(x_test);

figure(3)
plot(x_test, error_10, 'k-', 'LineWidth', 2); hold on;
title('Error: f(x) - p_{10}(x)');
xlabel('x'); ylabel('Error');
grid on;
text(x_center, max(error_10)*0.9, sprintf('Max Error = %.4e', max_error_10), ...
    'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
hold off

figure(4)
plot(x_test, error_16, 'k-', 'LineWidth', 2); hold on;
title('Error: f(x) - p_{16}(x)');
xlabel('x'); ylabel('Error');
grid on;
text(x_center, max(error_16)*0.9, sprintf('Max Error = %.4e', max_error_16), ...
    'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
hold off
