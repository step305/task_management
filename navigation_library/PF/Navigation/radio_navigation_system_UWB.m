function [fltr, robot_state] = radio_navigation_system_UWB(fltr, upsilon, dtheta, Ranges, Anchors, update, dt)
%   Алгоритм системы навигации колесного робота на базе локальной
%   радионавигационной системы UWB и фильтра частиц
%
%   Входные аргументы:
%   fltr - структура фильтра частиц
%   upsilon - управляющее воздействие робота - скорость продольного
%   перемещения
%   dtheta - управляющее воздействие робота - угловая скорость курсового
%   разворота
%   Ranges - измерения расстояний до базовых станций системы UWB
%   Anchors - координаты базовых станций системы UWB
%   update - индикатор доступности измерения от приемника UWB
%   dt - шаг интегрирования
%
%   Выходные аргументы:
%   fltr - структура фильтра частиц
%   robot_state - оценка вектора состояния робота

%% Предсказание/коррекция
fltr = pf_predict_update(fltr, upsilon, dtheta, Anchors, Ranges, update, dt);

%% Ресемплинг частиц
fltr = pf_resample(fltr, 0.75);

%% Состояние 
robot_state = pf_get_state(fltr, 2);

end



