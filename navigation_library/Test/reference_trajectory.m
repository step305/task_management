function [dThe_ref, dUps_ref, Rn_ref, Vb_ref, Euler_ref, Lat_ref, Lon_ref, dThet_ref, sample_rate] = reference_trajectory(sim_time)

% Параметры глобальной системы координат
a = 6378137;
e = 0.081819190842621;
lat_0 = 0.9;
lon_0 = 0.6;
alt_0 = 0.0;

% Траекторные параметры
s = 50;
fx = 0.05;
fy = 0.025;
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
Vb_ref = zeros(Nsamp,3);
Euler_ref = zeros(Nsamp,3);
dThe_ref = zeros(Nsamp,3);
dUps_ref = zeros(Nsamp,3);
Lat_ref = zeros(Nsamp,3);
Lon_ref = zeros(Nsamp,3);
dThet_ref = zeros(Nsamp,1);

% Начальные значения измерений датчиков угловой скорости и акселерометров
dthe = [0; 0; 0];
dups = [0; 0; 9.81e-4];

% Вектор ускорения свободного падение в навигационной системе координат
g = [0; 0; -9.81];

sample_cnt = 0;
t = 0;
for i=1:Nsim
    xr  = s*sin(fx*t*dt+fix);
    yr  = s*sin(fy*t*dt+fiy);
    dxr = fx*s*cos(fx*t*dt+fix);
    dyr = fy*s*cos(fy*t*dt+fiy);
    ur = sqrt(dxr^2+dyr^2);
    ddxr  = -fx^2*s*sin(fx*t*dt+fix);
    ddyr  = -fy^2*s*sin(fy*t*dt+fiy);
    dTr = (ddyr*dxr-ddxr*dyr)/(dxr^2+dyr^2);
    Tr = atan2(dyr,dxr);
    dur = -(s^2*(fx^3*sin(2*fix + 2*fx*t*dt) + ...
        fy^3*sin(2*fiy + 2*fy*t*dt)))/(2*(s^2*(fx^2*cos(fix + ...
        fx*t*dt)^2 + fy^2*cos(fiy + fy*t*dt)^2))^(1/2));
    fb = [dur;ur*dTr;0]+g;
    wb = [0;0;dTr];
    dthe = dthe+wb*dt;
    dups = dups+fb*dt;
    if ((mod(i,(sample_rate/dt)) == 0) || (t==0))
        sample_cnt = sample_cnt+1;
        Vb_ref(sample_cnt,:) = [ur;0;0];
        Rn_ref(sample_cnt,:) = [xr;yr;0];
        Euler_ref(sample_cnt,:)  = [Tr;0;0];
        dThe_ref(sample_cnt,:) = dthe;
        dUps_ref(sample_cnt,:) = dups;
        [Lat_ref(sample_cnt,:), Lon_ref(sample_cnt,:)] = ned_to_geo([xr; yr; 0], lat_0, lon_0, alt_0, a, e);
        dThet_ref(sample_cnt,:) = dTr;
        dthe = [0;0;0];
        dups = [0;0;0];
    end
    t = t+1;
end


