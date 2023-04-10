function p = make_ranging_lines (rb, xv, len)
% Вспомогательная функция для отображения измерений радара

if isempty(rb), p =[]; return, end
lnes(1,:) = zeros(1,len)+ xv(1);
lnes(2,:) = zeros(1,len)+ xv(2);
lnes(3,:) = zeros(1,len)+ xv(1)+rb(:,2)';
lnes(4,:) = zeros(1,len)+ xv(2)+rb(:,3)';
p = line_plot_conversion(lnes);
end
