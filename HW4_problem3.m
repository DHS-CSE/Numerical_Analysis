function PECE_adams_bashforth_moulton
    % Define parameters
    h = 1/16;
    t0 = 1;
    tf = 2;
    N = round((tf - t0)/h);
    t = t0:h:tf;
    y = zeros(1, N+1);
    y(1) = 1;  % Initial condition

    % Define RHS function
    f = @(t, y) 1/(y^2) - y*t;

    % Generate y(2), y(3), y(4) using RK4
    for i = 1:3
        k1 = f(t(i), y(i));
        k2 = f(t(i) + h/2, y(i) + h*k1/2);
        k3 = f(t(i) + h/2, y(i) + h*k2/2);
        k4 = f(t(i) + h, y(i) + h*k3);
        y(i+1) = y(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    % Apply PECE (Adams-Bashforth + Adams-Moulton)
    for i = 4:N
        % Get function values
        f_k   = f(t(i), y(i));
        f_k1  = f(t(i-1), y(i-1));
        f_k2  = f(t(i-2), y(i-2));
        f_k3  = f(t(i-3), y(i-3));

        % Predictor (Adams-Bashforth 4)
        y_pred = y(i) + h * (55*f_k - 59*f_k1 + 37*f_k2 - 9*f_k3)/24;

        % Evaluate f at predicted point
        f_kp1 = f(t(i+1), y_pred);

        % Corrector (Adams-Moulton 4)
        y(i+1) = y(i) + h * (9*f_kp1 + 19*f_k - 5*f_k1 + f_k2)/24;
    end

    % Final result
    fprintf("Approximate value of y(2) = %.8f\n", y(end));

    % Optional: plot
    plot(t, y, '-o');
    title('PECE Method (Adams-Bashforth-Moulton)');
    xlabel('t'); ylabel('y(t)');
    grid on;
end

PECE_adams_bashforth_moulton
