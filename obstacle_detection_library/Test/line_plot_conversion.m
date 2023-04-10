function p= line_plot_conversion (lne)
% Вспомогательная функция для отображения измерений радара

len= size(lne,2)*3 - 1;
p= zeros(2, len);
p(:,1:3:end)= lne(1:2,:);
p(:,2:3:end)= lne(3:4,:);
p(:,3:3:end)= NaN;
end