clc
clear
close all

%% Datos problema
E = 21e9;
t = 1;
v = 0.2;
rho = 2400;
%Plain stress
C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

W = 1;
L = 1;

bordes = [0 0;
          L 0;
          L W;
          0 W];

divx = 2;
divy = 2;

%% Mallador
elemtype = upper('Q8');
mesh = Q8(bordes,divx,divy,C);

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

%mesh.stress = mesh.funciones.stress(mesh,U);

% print_Stress_plane(mesh,U)


%% Funciones

function [plane] = plane_type(type)
    if strcmp(type, 'STRESS')
        plane.C = @(E,v,t) t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
        plane.Czz = @(E,v,t) zeros(1,3);
    elseif strcmp(type, 'STRAIN')
        plane.C = @(E,v,t) t*E/((1+v)*(1-2*v))*[1-v v 0;v 1-v 0;0 0 (1-2*v)/2];
        plane.Czz = @(E,v,t) E*v/(1+v)/(1-2*v)*[1,1,0];
    else
        warning('no se selecciono plane Stress o Strain')
    end
end

function [funciones] = elemento(elemType)
    switch elemType
        case 'Q4'
            libreria = [pwd '/Q4'];%asegurarse de cargar las funciones
            addpath(libreria);

            funciones.K = @K_Q4;
            funciones.mallador = @mallador_cuadrado_Q4;
            funciones.cargas = @cargas_Q4;
            funciones.stress = @stress_Q4;
        case 'Q8'
            libreria = [pwd '/Q8'];%asegurarse de cargar las funciones
            addpath(libreria);

            funciones.K = @K_Q8;
            funciones.mallador = @mallador_cuadrado_Q8;
            funciones.cargas = @cargas_Q8;
            funciones.stress = @stress_Q8;
         case 'CST'
            libreria = [pwd '/CST'];%asegurarse de cargar las funciones
            addpath(libreria);

            funciones.K = @K_CST;
            funciones.mallador = @mallador_triang_CST;
            funciones.cargas = @cargas_CST;
            funciones.stress = @stress_CST;
         case 'LST'
            libreria = [pwd '/LST'];%asegurarse de cargar las funciones
            addpath(libreria);

            funciones.K = @K_LST;
            funciones.mallador = @mallador_triang_LST;
            funciones.cargas = @cargas_LST;
            funciones.stress = @stress_LST;
        otherwise
            warning('elemento no valido')
    end
end

