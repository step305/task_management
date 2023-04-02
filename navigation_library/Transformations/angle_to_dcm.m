function Cbn = angle_to_dcm(rz, ry, rx)
%  Cbn = angle_dcm(rz, ry, rx)
%  Преобразует углы Эйлера в матрицу направляющих косинусов
%
%   Входные аргументы:
%   rz, ry, rx - углы Эйлера вокру осей Z, Y и X соответственно
%
%   Выходные аргументы:
%   Cbn -  матрица направляющих косинусос [3,3]

sz = sin(rz);
cz = cos(rz);
sy = sin(ry);
cy = cos(ry);
sx = sin(rx);
cx = cos(rx);

Cx = [1 0 0; 0 cx sx; 0 -sx cx];
Cy = [cy 0 -sy; 0 1 0; sy 0 cy];
Cz = [cz sz 0; -sz cz 0; 0 0 1];

Cbn = Cx * Cy * Cz;
