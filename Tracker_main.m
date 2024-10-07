%% Входные параметры
OMEGA_0 = deg2rad(15); % Долгота восходящего узла 
i = deg2rad(97.626); % Наклонение
u0 = deg2rad(0); % Начальный аргумент широты КА
H = 561.4 * 10^3; % Высота орбиты
C = 4; % Количество оборотов КА
N = 5000; % Общее количество рисуемых точек

%% Отрисовка

[phi_1, lambda_1] = coordinates_calculation(OMEGA_0, i, u0, H, C, N);

[phi_2, lambda_2] = coordinates_calculation(OMEGA_0, i, u0 + pi, H, C, N);

geoplot(phi_1, lambda_1, '.', phi_2, lambda_2, '.r');
