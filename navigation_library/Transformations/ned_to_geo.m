function [lat, lon, alt] = ned_to_geo(xn, lat0, lon0, alt0, a, e)
%  [lat, lon, alt] = ned2geo(xn, lat0, lon0, alt0, a, e)
%  Преобразует вектор координат из локальной системы координат NED в географическую
%  систему координат
%
%   Входные аргументы:
%   xn -  Координаты точки в системе NED
%   lat0, lon0, alt0 - Широта, долгота и высота начальной точки локальной
%   системы координат
%   a, e  -  Параметры эллипсоида географической системы координат
%
%   Выходные аргументы:
%   lat, lon, alt - Широта, долгота и высота точки в географической системе
%   координат 

xe0  = geo_to_ecef(lat0, lon0, alt0, a, e);
Cen = [...
    -sin(lat0) * cos(lon0) -sin(lat0) * sin(lon0)  cos(lat0)
                -sin(lon0)              cos(lon0)          0
    -cos(lat0) * cos(lon0) -cos(lat0) * sin(lon0) -sin(lat0)
    ];
xe = xe0 + Cen' * xn;
[lat, lon, alt] = ecef_to_geo(xe, a, e);
end