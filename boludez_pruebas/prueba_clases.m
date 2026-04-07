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

%% Q4
cant_puntos = 4;

% creo el mesh y calculo la K
mesh = Q4(bordes,divx,divy,C);

% hacer la carga
carga = zeros(4,4);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];

llamar_funciones(mesh,carga,volumen,50);

%% Q8
cant_puntos = 8;

% creo el mesh y calculo la K
mesh = Q8(bordes,divx,divy,C);

% hacer la carga
carga = zeros(4,6);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];

llamar_funciones(mesh,carga,volumen,50);

%% CST
cant_puntos = 3;

% creo el mesh y calculo la K
mesh = CST(bordes,divx,divy,C);

% hacer la carga
carga = zeros(3,4);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];

llamar_funciones(mesh,carga,volumen,50);

%% LST
cant_puntos = 6;

% creo el mesh y calculo la K
mesh = LST(bordes,divx,divy,C);

% hacer la carga
carga = zeros(3,6);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];

llamar_funciones(mesh,carga,volumen,50);

%% funciones
function [] = llamar_funciones(mesh,carga,volumen,porcentaje)
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
end