function move_frame(Cbn, Rn, handle1, handle2, handle3, a1, a2, a3)
% Вспомогательная функция для анимации перемещения препятствия

a1_ = Cbn*a1+Rn;
a2_ = Cbn*a2+Rn;
a3_ = Cbn*a3+Rn;
s_ = Rn;
set(handle1,'Xdata',[s_(1) a1_(1)]);
set(handle1,'Ydata',[s_(2) a1_(2)]);
set(handle1,'Zdata',[s_(3) a1_(3)]);
set(handle2,'Xdata',[s_(1) a2_(1)]);
set(handle2,'Ydata',[s_(2) a2_(2)]);
set(handle2,'Zdata',[s_(3) a2_(3)]);
set(handle3,'Xdata',[s_(1) a3_(1)]);
set(handle3,'Ydata',[s_(2) a3_(2)]);
set(handle3,'Zdata',[s_(3) a3_(3)]);
end