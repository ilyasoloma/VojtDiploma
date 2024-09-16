R_earth = 6371 * 10^3; %Радиус Земли
mu = 3.986*10^(14); % гравитационный параметр Земли, м^3/с^2
OMEGA_0 = 15 * pi / 180; 
i = 97.626 * pi / 180;
% tau = 0;
omega_earth = 0.71921 * 10^(-4);

N = 2000;

H = 561.4 * 10^3;


r = R_earth + H; %Радиус апогея

a = r;

T_star = 2 * pi * sqrt(a^3 / mu);

p = a;

d_OMEGA = (-35.062 / 60) * (R_earth / p)^2;

dt = 7 * T_star / N;

OMEGA = zeros(1, N);
phi = zeros(1, N);
argument_phi = zeros(1, N);
argument_sin = zeros(1, N);
phi_deriv = zeros(1, N);
lambda = zeros(1, N);

for j = 1:N
    
    t = j * dt;
    
    OMEGA(j) = OMEGA_0 + (t / T_star) * d_OMEGA;
    
    % Широта
    argument_sin(j) = 2*pi*t/T_star;
    argument_phi(j) = sin(i) * sin(argument_sin(j));
    phi(j) = asin(argument_phi(j));
    
    phi_deriv(j) = (sin(i) * cos(argument_sin(j)) * 2 * pi / T_star) / (sqrt(1 - argument_phi(j)^2));
    
    % Долгота
    lambda_check = OMEGA(j) + atan(cos(i) * tan(2*pi*t/T_star)) - omega_earth * t + d_OMEGA / T_star;
    
    if phi_deriv(j) < 0
        lambda_check = lambda_check + pi;
    elseif phi_deriv(j) > 0 && phi(j) < 0
        lambda_check = lambda_check + 2*pi;
    end
    
    while lambda_check < 0 || lambda_check > 2 * pi
        if lambda_check < 0
            lambda_check = lambda_check + 2 * pi;
        elseif lambda_check > 2 * pi
            lambda_check = lambda_check - 2 * pi;
        end
    end
    
    lambda(j) = lambda_check;
end

lambda_debug = lambda ./ pi;

phi = rad2deg(phi);

lambda = rad2deg(lambda);

geoplot(phi, lambda, '.');
