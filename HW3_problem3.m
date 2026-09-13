% 3D composite trapezoid rule

syms xs ys zs
f_sym = xs * cos(xs*ys) + zs.^2;
exact_int = int(int(int(f_sym, zs, 0, pi), ys, 0, pi), xs, 0, pi);
disp('Exact_int =')
disp(double(exact_int))

a = 0; b = pi;

N = 10;
tol = 1e-2;

while true
    h = (b - a) / N;
    x = linspace(a, b, N+1);
    y = linspace(a, b, N+1);
    z = linspace(a, b, N+1);
    dx = h;
    dy = h;
    dz = h;

    wx = ones(1, N+1); wx([1 end]) = 0.5;
    wy = ones(1, N+1); wy([1 end]) = 0.5;
    wz = ones(1, N+1); wz([1 end]) = 0.5;

    I = 0;

    for k = 1:N+1
        zk = z(k);
        Iz = 0;
        
        for j = 1:N+1
            yj = y(j);
            Iy = 0;
            
            for i = 1:N+1
                xi = x(i);
                f_vals = xi * cos(xi * yj) + zk^2;
                Iy = Iy + f_vals * wx(i);
            end
            
            Iy = Iy * dx;
            Iz = Iz + Iy * wy(j);
        end
        
        Iz = Iz * dy;
        I = I + Iz * wz(k);
    end
    
    I = I * dz;
    
    d2f_dx2 = diff(f_sym, xs, 2);
    d2f_dy2 = diff(f_sym, ys, 2);
    d2f_dz2 = diff(f_sym, zs, 2);
    laplacian_f = d2f_dx2 + d2f_dy2 + d2f_dz2;
    
    [Xg, Yg, Zg] = ndgrid(x, y, z);
    lap_fn = matlabFunction(laplacian_f, 'Vars', [xs, ys, zs]);
    lap_vals = lap_fn(Xg, Yg, Zg);
    f_dd_max = max(abs(lap_vals));

    error = (b - a)^3 * h^2 / 12 * f_dd_max;
    if error > tol
        N = N+1;
    else
        approximation = I;
        break;
    end
end

fprintf('Approximate Integral: %.6f\n', approximation);
fprintf('Final step size h = %.6f with N = %d\n', h, N);
%% monte carlo method

syms x y z
f = x * cos(x*y) + z.^2;
exact_int = int(int(int(f, z, 0, pi), y, 0, pi), x, 0, pi);
disp('Exact_int =')
disp(double(exact_int))

f = matlabFunction(f);
a = 0;
b = pi;
iteration = 1;
N_total = zeros(1, iteration);
k = 1;
tol = 1e-2;

while k <= iteration
    N = 10000000;
    while true
        temp_xyz = (b - a) * rand(N, 3);
        values = zeros(N, 1);
        var = zeros(N, 1);
        for i = 1:N
            values(i) = f(temp_xyz(i, 1), temp_xyz(i, 2), temp_xyz(i, 3));
        end
        monte_integration = mean(values) * (b-a).^3;
        for j = 1:N
            var(j) = (values(j) - mean(values)).^2;
        end
        standard_deviation = sqrt(mean(var));
        err = standard_deviation / sqrt(N);
    
        if err > tol
            N = N + 1;
        else
            break;
        end
    end
    N_total(k) = N;
    k = k + 1;
end
fprintf('Approximate Monte Integral: %.6f\n', monte_integration);
fprintf('Final nodes number N = %d\n', mean(N_total));

