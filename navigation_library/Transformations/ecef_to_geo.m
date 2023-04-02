function [lat, lon, alt] = ecef_to_geo(xe, a, e)
%  [lat, lon, alt] = ecef2geo(xe, a, e)
%  Преобразует координаты точки их глобальной декартовой системы координат ECEF
%  в глобальную географическую систему координат (СК)
%
%   Входные аргументы:
%   xe -  Координаты  точки в СК ECEF
%   a, e  -  Параметры эллипсоида географической СК
%
%   Выходные аргументы:
%   lat, lon, alt - Широта, долгота и высота точки в географической СК

lon = atan2(xe(2, 1), xe(1, 1));
alt = 0;
Re = a;
p = sqrt(xe(1, 1) ^ 2 + xe(2, 1) ^ 2);

for i=1:25
    sin_lat = xe(3, 1)/((1 - e ^ 2) * Re + alt);
    lat = atan((xe(3, 1) + e ^ 2 * Re * sin_lat) / p);
    Re = a / sqrt(1 - e ^ 2 * sin(lat) ^ 2);
    alt = p / cos(lat) - Re;
end
end