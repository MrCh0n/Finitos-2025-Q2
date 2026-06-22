clc
clear
close all

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos
L = 30e-3;%m
W = L;%m
r = 10e-3;

E = 10e9;%Gpa
v = 0.3;
t = 1;

F = 30e3;% Carga total [N]

%Stress
C = E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
Czz = zeros(1,3);

bordes = [0 0;
          L 0;
          L W;
          0 W];

divx = 4;
divy = 2;
%% Mesh
[nodos, elems, bordes] = mallador_cuadrado_Q4(bordes, divx, divy);

nnod = size(nodos,1);
ndof = nnod*2;
nelem = size(elems,1);

dofs = reshape(1:ndof,2,[])';

free = true(nnod, 2);
%empotrado a izquierda, abajo, derecha
empotrado = [bordes.lado_41, bordes.lado_12, bordes.lado_23];
free(empotrado,:) = false;

free = reshape(free', 1, []);
%% Crear K y R
dofselem = 8;
crearK = @crearK_Q4;

K = armar_K_2D(nodos,elems,dofs,dofselem,C,crearK);
%% Hacer el aumento de lagrange
parte_superior = dofs(bordes.lado_34,2);%y de los nodos superiores

nlag = size(parte_superior,1);
free_lag = true(1,nlag);
G = zeros(ndof,nlag);
profundidad = 8e-3;
for i = 1:nlag
    G(parte_superior(i), i) = -1;
    abs(nodos(bordes.lado_34(i),1)-L/2)
    if abs(nodos(bordes.lado_34(i),1)-L/2)>profundidad
        free_lag(i) = false;
    end
end
%% Crear R y g
R = zeros(ndof,1);
R(parte_superior([ceil(end/2) ceil(end/2)-1 ceil(end/2)+1])) = -F/3;
g = zeros(nlag,1);
g(free_lag) = profundidad-(nodos(bordes.lado_34(free_lag),1)-L/2).^2/(2*r);
%% Calculo
frees = [free, free_lag];

dofs_totales = ndof+nlag;
K_aum = [K G; G' zeros(nlag,nlag)];
R_aum = [R;g];

Kr = K_aum(frees, frees);
Rr = R_aum(frees);

U = zeros(dofs_totales,1);
U(frees) = Kr\Rr;

u = U(1:ndof);
lambda = U(ndof+1:end);