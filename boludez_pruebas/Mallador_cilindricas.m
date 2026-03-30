clc
clear

addpath(genpath(pwd+"/../Libreria_elementos"))
a = 1;
b = 1;
bordes = [0 0;
          a 0;
          a b;
          0 b];
divx = 2;
divy = 2;
[nodos, elems] = mallador_cuadrado_Q4(bordes,divx,divy);

nnod = size(nodos,1);

R = 1;
titamax=40;
zmin = R*cosd(titamax);
for i = 1:nnod
    x = nodos(i,1);
    y = nodos(i,2);

    tita = titamax*(x/a);

    z = R*cosd(tita) - zmin;
    aux_nodos(i,:) = [x y z];
end