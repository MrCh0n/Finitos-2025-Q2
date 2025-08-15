clc 
clear

mesh = struct('coord',[],'elem',[],'prop',[],'BC',[],'Fnod',[],'Q', []);

E = 70e9;
A = 500e-6;

mesh.coord = [0 0;
              10 0;
              20 0;
              30 0;
              10 10;
              20 10]/10*3.048;
angulo = 30;

Q = [cosd(angulo) sind(angulo); -sind(angulo) cosd(angulo)];
mesh.Q = Q;

for i = 1:size(mesh.coord,1)
    mesh.coord(i,:) = Q*mesh.coord(i,:)';
end

mesh.elem = [1 2;
             2 3;
             3 4;
             4 6;
             5 6;
             1 5;
             3 6;
             2 6;
             2 5];

for i = 1:size(mesh.elem,1)
    mesh.prop(i,:) = [E A];
end

mesh.BC = [1 2 8];

P1 = Q*[0 -400]';
P2 = Q*[0 -300]';
P3 = Q*[0 -100]';

mesh.Fnod = [1 P1(1);
             2 P1(2);
             3 P2(1);
             4 P2(2);
             5 P3(1);
             6 P3(2)];

mesh.Fnod(:,2) = mesh.Fnod(:,2)/0.225;
barras(mesh);