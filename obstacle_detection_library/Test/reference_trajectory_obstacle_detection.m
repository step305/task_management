function [Rn_ref, Euler_ref, sample_rate] = reference_trajectory_obstacle_detection(sim_time)
%   Входные аргументы:
%   sim_time - время моделирования, секунды
%
%   Выходные аргументы:
%   Rn_ref -  координаты платформы
%   Euler_ref -  углы ориентации платформы
%   sample_rate -  шаг моделирования 

% Траекторные параметры
s = 30;
fx = 0.1;
fy = 0.05;
fix = 0;
fiy = 0;

% Шаг моделирование и шаг модельного опроса датчиков
dt = 1/10000;
sample_rate = 1/100;

% Время моделирования
Nsamp = sim_time*(sample_rate/dt);

% Логи
Nsim = sim_time/dt;
Rn_ref = zeros(Nsamp,3);
Euler_ref = zeros(Nsamp,3);

sample_cnt = 0;
t = 0;
for i=1:Nsim
    xr  = s*sin(fx*t*dt+fix);
    yr  = s*sin(fy*t*dt+fiy);
    dxr = fx*s*cos(fx*t*dt+fix);
    dyr = fy*s*cos(fy*t*dt+fiy);
    Tr = atan2(dyr,dxr);
    if ((mod(i,(sample_rate/dt)) == 0) || (t==0))
        sample_cnt = sample_cnt+1;
        Rn_ref(sample_cnt,:) = [xr;yr;0];
        Euler_ref(sample_cnt,:)  = [Tr;0;0];
    end
    t = t+1;
end


