R_earth = 6371 * 10^3; %Радиус Земли
mu = 3.986*10^(14); % гравитационный параметр Земли, м^3/с^2
OMEGA_0 = 15 * pi / 180; 
i = 97 * pi / 180;
% tau = 0;
omega_earth = 0.71921 * 10^(-4);

N = 500;

H = 561.4 * 10^3;


r = R_earth + H; %Радиус апогея

a = r;

T_star = 2 * pi * sqrt(a^3 / mu);

p = a;

d_OMEGA = (-35.062 / 60) * (R_earth / p)^2;



% d_omega = (-17.525 / 60) * (1 - 5 * (cos(i))^2) * (R_earth / p)^2;

% omega = omega_0 + (t / T_star) * d_omega;

dt = T_star / N;

OMEGA = zeros(1, N);
phi = zeros(1, N);
lambda = zeros(1, N);

for j = 1:N
    
    t = j * dt;
    
    OMEGA(j) = OMEGA_0 + (t / T_star) * d_OMEGA;
    
    % Широта
    phi(j) = asin(sin(i) * sin(2*pi*t/T_star));
    % Долгота
    lambda(j) = OMEGA(j) + atan(cos(i) * tan(2*pi*t/T_star)) - omega_earth * t + d_OMEGA / T_star;
end

phi = rad2deg(phi);

lambda = rad2deg(lambda);

geoplot(phi, lambda);
