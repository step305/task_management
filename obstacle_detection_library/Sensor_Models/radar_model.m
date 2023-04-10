function [target_pos, range, azimuth, visible] = radar_model(rn_platform, Cbn_platform, rn_target, scan_angle, range_max, range_min)
% Базовая модель измерений радара или сканирующего дальномера
%
%   Входные аргументы:
%   rn_platform - Координаты платформы в навигационной СК
%   Cbn_platform - Матрицы направляющих косинусов ориентации платформы
%   rn_target  - Координаты препятствия в навигационной СК
%   scan_angle - Угол раствора сканирования радара или сканирующего
%   дальномера
%   range_max  - Максимальная рабочая дальность радара или сканирующего дальномера
%   range_min  - Минимальная рабочая дальность радара или сканирующего дальномера
%
%   Выходные аргументы:
%   target_pos - Зашумленные координаты препятствия 
%   range    - Дальность до препятствия
%   azimuth  - Азимутальный угол пеленга преяптствия в связанной с
%   платформой СК
%   visible  - Индикатор видимости препятствия

fb = Cbn_platform' * (rn_target - rn_platform);
azimuth = atan2(fb(2), fb(1));
range = norm(rn_target - rn_platform);

visible = false;

if (abs(azimuth) < scan_angle)
    visible = true;
end

if ((range > range_max) && (range < range_min))
    visible = false;
end

% Ошибка измерения угла
angle_noise = 0.05;
azimuth = azimuth + randn * angle_noise;

% Ошибка измерения дальности
range_noise = 0.2;
range = range + randn * range_noise;

target_pos = rn_platform + angle_to_dcm(azimuth, 0, 0)' * Cbn_platform * [range; 0; 0];
target_pos = target_pos(1:2, 1);

end
