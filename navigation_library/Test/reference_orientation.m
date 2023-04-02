function [dWb, q_ref, dt, q_init, Fb, Mb, fn, mn, bw] = reference_orientation
% [dWb, q, dt, q_init, fb, mb, fn, mn, bw] = reference_data_spin_cone
% Синтезирует тестовые данные для проверки функционирование алгоритмов
% определения ориентации 
% Опорное движение твердого тела задается с помощью алгоритма SPIN-CONE
% описывающего тело, вращающееся с постоянной угловой скоростью вокруг
% некоторой оси, которая в свою очередь описывает конус с постоянной
% скоростью прецесии
%
%   Выходные аргументы:
%   dWb -  Измерения датчиков угловой скорости (массив)
%   q  - Опорный кватернион ориентации (массив)
%   dt - Длина шага интегрирования 
%   q_init - Начальное значение опорного кватерниона ориентации
%   fb  - Измерения акселерометров (массив)
%   mb  - Измерения магнитометров (массив)
%   fn - Вектор ускорения свободного падение в навигационной системе
%   координат 
%   mn - Вектор магнитного поля в навигационной системе координат
%   bw - Опорные сдвиги нулей датчиков угловой скорости

%%
clear;
clc;

%%
rng('shuffle');

%% Параметры моделирования
% Время моделирования, сек.
tsim = 360;
% Шаг интегрирования, сек.
dt = 1e-2;

%% Кинематические параметры
% Угловая скорость прецесии, рад/сек.
wc = 1e-1;
% Угловая скорость вращения, рад/сек.
ws = 1e0;
% Угол между осями вращения и прецесии, рад.
beta = pi/4;

%% Шумы датчиков
nw = 1e-8; % датчики угловой скорости 
nf = 1e-6; % акселерометры
nm = 1e-6; % магнитометры

%% Сдвиги нолей датчиков
bdw = [1e-4; -1e-4; 2e-4]; % датчики угловой скорости 
bfb = [0e0; 0e0; 0e0];     % акселерометры
bmb = [0e0; 0e0; 0e0];     % магнитометры

%% Шумы случайного блуждания сдвигов нолей датчиков угловой скорости 
rwn = 1e-10; 

%% Опорные векторы в навигационной системе координат 
fn = [0;0;-1]; % вектор ускорения свободного падения
mn = [1;0;0];  % вектор магнитного поля

%%
dWb = zeros(tsim/dt,3);
q_ref   = zeros(tsim/dt,4);
angles = zeros(tsim/dt,3);
Fb = zeros(tsim/dt,3);
Mb = zeros(tsim/dt,3);
bw = zeros(tsim/dt,3);

%% Моделирование
time = 0;
cnt = 0;
[iWb_, Cbn] = spin_cone( wc, ws, beta, time );
q_init = dcm_to_quat(Cbn);
while (time < tsim)
    time = time+dt;
    cnt = cnt+1;
    bdw = bdw+randn(3,1)*rwn;
    [ iWb, Cbn, rz, ry, rx ] = spin_cone( wc, ws, beta, time );
    da = iWb-iWb_;
    dWb(cnt,:) = da+randn(3,1)*sqrt(nw)+bdw;
    q_ref(cnt,:) = dcm_to_quat(Cbn);
    angles(cnt,:) = [rz, ry, rx];
    Fb(cnt,:) = Cbn'*fn+randn(3,1)*sqrt(nf)+bfb;
    Mb(cnt,:) = Cbn'*mn+randn(3,1)*sqrt(nm)+bmb;
    bw(cnt,:) = bdw;
    iWb_ = iWb;
end
end


function [iWb, Cbn, psi, theta, phi] = spin_cone(wc, ws, beta, t)

% Углы ориентации
phi = (ws-wc*cos(beta))*t;
theta = pi/2-beta;
psi = -wc*t;

% Матрица направляющих косинусов 
Cbn = zeros(3);
Cbn(1,1) = cos(theta)*cos(psi);
Cbn(1,2) = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
Cbn(1,3) = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
Cbn(2,1) = cos(theta)*sin(psi);
Cbn(2,2) = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
Cbn(2,3) = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
Cbn(3,1) = -sin(theta);
Cbn(3,2) = sin(phi)*cos(theta);
Cbn(3,3) = cos(phi)*cos(theta);

% Интеграл угловой скорости
iWb = zeros(3,1);
if norm(ws) > 0
    iWb(1,1) = ws*t;
    iWb(2,1) = ((wc*sin(beta))/(ws-wc*cos(beta)))*cos(phi);
    iWb(3,1) = -((wc*sin(beta))/(ws-wc*cos(beta)))*sin(phi);
end
end


