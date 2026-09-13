L = 0.5;
rho = 1.0;
x0 = 0.5;

g = @(x) cos(x);

max_iter = 100;
tolerance = 1e-6;

x_k = x0;
x_star = g(x0);

x_vals = zeros(1, max_iter);
errors = zeros(1, max_iter);
rate_of_convergence = zeros(1, max_iter);
x_vals(1) = x_k;

for k = 1:max_iter
    x_k1 = g(x_k);
    x_vals(k+1) = x_k1;
    errors(k+1) = abs(x_k1 - x_star);

    if k > 1
        rate_of_convergence(k) = errors(k+1) / errors(k);
    end

    % 수렴 조건 체크
    if abs(x_k1 - x_k) < tolerance
        break;
    end
    x_k = x_k1;
end

% 결과 출력
disp('고정점 x* = ');
disp(x_k1);
disp('반복 횟수: ');
disp(k);

% 수렴 속도 확인 (그래프)
disp('수렴 속도 L (e_{k+1} / e_k):');
disp(rate_of_convergence(2:k));  % 첫 번째 값 제외

% 수렴 속도 그래프
figure;
plot(2:k, rate_of_convergence(2:k), 'o-');
xlabel('Iteration');
ylabel('Rate of Convergence e_{k+1} / e_k');
title('수렴 속도 측정');
grid on;