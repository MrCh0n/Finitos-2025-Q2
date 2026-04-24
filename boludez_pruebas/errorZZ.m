clc
clear

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos problema
E = 200e9;
t = 1;
v = 0.3;
%Plain stress
C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

a = 1;

bordes = [0 0;
          a 0;
          a a;
          0 a];

divx = 10;
divy = 10;

%% mashado
cant_puntos = 4;

% creo el mesh y calculo la K
mesh = Q4(bordes,divx,divy,C);

% hacer la carga
carga = zeros(4,4);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];

mesh.armar_K;

% i puede ser el elemento justo que se quiere cargar, puse todos porque
% quiero hacer una fuerza volumetrica solo
for i = 1:mesh.counts.nelem
   mesh.armar_R(i,carga,volumen,true);
end

% Condiciones de borde
% talvez se puede hacer mas facil, es empotrados en la parte de abajo
borde_inferior = mesh.bordes.lado_12;
mesh.cond_borde(borde_inferior,3);

% Calculo y plot
mesh.calc_U;

%Dibuja ~50% de escala del modelo
%mult = mesh.dibujar(porcentaje);
hold off

mesh.calc_errorzz;
mesh.error.zz
% dofselem = 8;
% suavizado = zeros(mesh.counts.nnod,3);%exx,eyy,exy
% deformaciones = zeros(mesh.counts.nnod,3*cant_puntos);
% cont = zeros(mesh.counts.nnod,3);
% 
% for i = 1:mesh.counts.nelem
%     nodosid = mesh.elems(i,:);
% 
%     coord = mesh.nodos.coordenadas(nodosid,:);
% 
%     dir = reshape(mesh.nodos.dofs(nodosid,:)',1,[]);
% 
%     Uel = mesh.U(dir);
% 
%     epsilon_el = deformaciones_Q4(coord, Uel);
% 
%     deformaciones(i,:) = epsilon_el';
% 
%     suavizado(nodosid,:) = suavizado(nodosid,:) + reshape(epsilon_el,3,[])';
%     cont(nodosid,:) = cont(nodosid,:) + 1;
% end
% 
% suavizado = suavizado./cont;
% 
% mesh.error.U = 0;
% mesh.error.E = 0;
% for i = 1:mesh.counts.nelem
%     nodosid = mesh.elems(i,:);
% 
%     coord = mesh.nodos.coordenadas(nodosid,:);
% 
%     e_el = deformaciones(i,:);
%     e2_el = reshape(suavizado(nodosid,:)',1,[])-e_el;
% 
%     Uel = energia_Q4(coord, e_el, mesh.material.C);
% 
%     Eel = energia_Q4(coord, e2_el, mesh.material.C);
% 
%     mesh.error.U = mesh.error.U + Uel;
%     mesh.error.E = mesh.error.E + Eel;
% end
% mesh.error.zz = sqrt(mesh.error.E/(mesh.error.E+mesh.error.U));
