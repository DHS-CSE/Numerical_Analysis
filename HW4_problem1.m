% Approximation with Forward Euler Method
h = 1/16;
num_points = int64(1/h) + 1;
t = linspace(0, 1, num_points);
N = length(t);
y1_euler = zeros(1, N);
y2_euler = zeros(1, N);
y1_euler(N) = 2;
y2_euler(N) = 1;

for i=N:-1:2
    y1_euler(i-1) = y1_euler(i) - h*y2_euler(i);
    y2_euler(i-1) = y2_euler(i) - h*((1-t(i))*y1_euler(i)+1)/(1+t(i)).^2;
end

% Backward Euler Method
y1_back = zeros(1, N);
y2_back = zeros(1, N);
y1_back(N) = 2;
y2_back(N) = 1;

for i=N:-1:2
    y2_old = y2_back(i);
    while true
        y1_prev = y1_back(i) - h * y2_old;
        f_val = ((1 - t(i-1)) * y1_prev + 1) / (1 + t(i-1))^2;
        y2_new = y2_back(i) - h * f_val;
    
        if abs(y2_new - y2_old) < 1e-12
            y2_back(i-1) = y2_new;
            y1_back(i-1) = y1_back(i) - h*y2_back(i-1);
            break;
        end  
        y2_old = y2_new;
    end
end

figure
plot(t, y1_euler, 'b-o', 'DisplayName', 'Euler Method');
hold on;
plot(t, y1_back, 'r-*', 'DisplayName', 'Backward Euler Method');
xlabel('t');
ylabel('y(t)');
legend;
title('Solution of terminal value problem by Euler and Backward Euler methods');
grid on;