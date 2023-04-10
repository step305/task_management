function x = pf_get_state(fltr, type)
% Оценка вектора состояния фильтра частциц

switch type
    case 1 % состояние частицы с наибольшим весом
        idx = fltr.w == max(fltr.w);
        x = mean(fltr.p(:,idx),2);
    case 2 % Среднее состояние по всем частицам
        wp = fltr.p .* repmat(fltr.w, size(fltr.p, 1), 1);
        x = sum(wp, 2);
        x(3, 1) = pi2pi(x(3, 1));
end

end