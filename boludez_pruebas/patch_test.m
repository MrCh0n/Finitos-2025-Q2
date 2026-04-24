clc
clear
close all

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos problema
E = 21e9;
t = 1;
v = 0.2;
rho = 2400;

W = 1;
L = 1;

bordes = [0 0;
          L 0;
          L W;
          0 W];

divx = 2;
divy = 2;

%% Mallador
elemtype = upper('Q4');
mesh = Q4(bordes,divx,divy,E,v,t,"stress");

mesh.nodos.coordenadas(ceil(mesh.counts.nnod/2),:) = 1.1*mesh.nodos.coordenadas(ceil(mesh.counts.nnod/2),:);

switch elemtype
    case 'CST'
        F1 = 500e3;
        F2 = 1000e3;
    case 'Q4'
         F1 = 500e3;
        F2 = 1000e3;
    case 'LST'
        F1 = 1000e3;
        F2 = 4000e3;
    case 'Q8'
        F1 = 1000e3;
        F2 = 4000e3;
end



%% Condiciones de Borde
borde_izquierdo = mesh.bordes.lado_41;
mesh.cond_borde(borde_izquierdo,1); %fijo x lado izquierdo
mesh.cond_borde(borde_izquierdo(1),2);%el 1 va empotrado


%% Carga

R = zeros(mesh.counts.ndof,1);

%Borde derecho
borde_derecho = mesh.bordes.lado_23;
dof_borde = mesh.nodos.dofs(borde_derecho,1);%solo x

R(dof_borde(2:end-1)) = F2; 
R(dof_borde(1)) = F1; 
R(dof_borde(end)) = F1;

if strcmp(elemtype, 'Q8')
    R(dof_borde(3)) =  R(dof_borde(3))/2;
end

mesh.R = R;

%% Armar matriz

mesh.armar_K;

%% Calculos
mesh.calc_U;

%% Graficos Deformacion
mesh.dibujar(100);

%% Tensiones
mesh.calc_stress;
%mesh.stress = mesh.funciones.stress(mesh,U);

% print_Stress_plane(mesh,U)