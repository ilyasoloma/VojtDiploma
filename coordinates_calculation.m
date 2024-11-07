function [phi, lambda] = coordinates_calculation(OMEGA_0, i, u0, H, C, N)
%% Постоянные
R_earth = 6371 * 10^3; %Радиус Земли
mu = 3.986*10^(14); % гравитационный параметр Земли, м^3/с^2
omega_earth = 0.71921 * 10^(-4); % Угловая скорость вращения Земли

%% Расчет

r = R_earth + H; %Радиус орбиты

a = r; % Большая полуось

T_star = 2 * pi * sqrt(a^3 / mu); % Период обращения КА

p = a; % Фокальный параметр

d_OMEGA = (-35.062 / 60) * (R_earth / p)^2; % Прецессия линии узлов за период обращения КА



dt = C * T_star / N; % Временной промежуток между отрисованными точками

OMEGA = zeros(1, N); % Текущая долгота восходящего узла
phi = zeros(1, N); % Текущая широта
argument_phi = zeros(1, N); 
argument_sin = zeros(1, N);
phi_deriv = zeros(1, N);
lambda = zeros(1, N); % Текущая долгота

for j = 1:N
    
    t = j * dt; 
    
    u = ((2 * pi * t) / T_star) + u0;
    
    OMEGA(j) = OMEGA_0 + (t / T_star) * d_OMEGA;
    
    % Широта
    argument_sin(j) = u;
    argument_phi(j) = sin(i) * sin(argument_sin(j));
    phi(j) = asin(argument_phi(j));
    
    phi_deriv(j) = (sin(i) * cos(argument_sin(j)) * 2 * pi / T_star) / (sqrt(1 - argument_phi(j)^2));
    
    % Долгота
    lambda_check = OMEGA_0 + atan(cos(i) * tan(u)) - omega_earth * t + d_OMEGA / T_star;
%     lambda_check = asin((sin(OMEGA(j))*cos(u) + cos(OMEGA(j))*cos(i)*sin(u)) / (cos(phi(j))));

    if phi_deriv(j) < 0
        lambda_check = lambda_check + pi;
    elseif phi_deriv(j) > 0 && phi(j) < 0
        lambda_check = lambda_check + 2*pi;
    end
    
    lambda_check = mod(lambda_check, 2 * pi);
    
%     while lambda_check < 0 || lambda_check > 2 * pi
%         if lambda_check < 0
%             lambda_check = lambda_check + 2 * pi;
%         elseif lambda_check > 2 * pi
%             lambda_check = lambda_check - 2 * pi;
%         end
%     end
    
    lambda(j) = lambda_check;
end

phi = rad2deg(phi);

lambda = rad2deg(lambda);
end
