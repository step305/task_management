function detector = basic_obstacle_detector(detector, rn_platform, Cbn_platform, targets_range, targets_azimuth, targets_visible)
% Базовая логика детектора подвижных и неподвижных препятствий -
% подразумевает известными идентификаторы для каждого из препятствий
%
%   Входные аргументы:
%   detector - структура детектора препятствий с трекерами на базе EKF и IMM 
%   rn_platform -  координаты платформы
%   Cbn_platform - углы ориентации платформы
%   targets_range - дальности до препятствий
%   targets_azimuth - углы пеленга препятствий
%   targets_visible - индикаторы видимости препятствий
%
%   Выходные аргументы:
%   detector - структура детектора препятствий с трекерами на базе EKF и IMM 

ntargets = size(targets_visible, 1);

latch_hits = 3;
max_misses = 10;

for j=1:ntargets
    if strcmp(detector{j}.type, 'imm')
        detector{j} = imm_tracker(detector{j}, rn_platform, Cbn_platform, targets_range(j, :), targets_azimuth(j, :), targets_visible(j, 1));
    elseif strcmp(detector{j}.type, 'ekf')
        detector{j} = ekf_tracker(detector{j}, rn_platform, Cbn_platform, targets_range(j, :), targets_azimuth(j, :), targets_visible(j, 1));
    end

    if targets_visible(j, 1) == 1
        detector{j}.hits = detector{j}.hits + 1;
        detector{j}.misses = 0;
        if (detector{j}.hits) > latch_hits
            detector{j}.latched = 1;
        end
    else
        detector{j}.hits = 0;
        if detector{j}.misses < max_misses
            detector{j}.misses = detector{j}.misses + 1;
        else           
            if strcmp(detector{j}.type, 'imm')
                detector{j}.P  = diag([100, 100, 1, 1]);
                detector{j}.P1 = diag([100, 100, 1, 1]);
                detector{j}.P2 = diag([100, 100]);
                detector{j}.mu = [0.5, 0.5];
            elseif strcmp(detector{j}.type, 'ekf')
                detector{j}.P  = diag([100, 100, 1, 1]);
            end
            detector{j}.latched = 0;
        end
    end

end

end