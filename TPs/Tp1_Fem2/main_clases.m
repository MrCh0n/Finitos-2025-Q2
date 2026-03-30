% Main del trabajo practico 1

clc
clear all
close all

%% Datos del problema
E = 4.32e8; %psi
v = 0;

R = 25; %in
L = 50; %in
t = 0.25 ; %in
tita = 40; %grados

q = -90; %psi

%% Control
type = 2; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 16; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        mesh = Mindlin(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5 6];
        sym_xz = [2 4 6];
    case 2
        mesh = Degenerado(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5];
        sym_xz = [2 4];
end
%% Malla
%ya esta mallado necesito nodos para las BC
nodos = mesh.nodes;

%% Matriz de rigidez
mesh.armar_K;

%% Condiciones de borde
borde_AC = abs(nodos(:,1)) <1e-6; %x = 0
borde_CD = abs(nodos(:,2)) <1e-6; %y = 0
borde_AB = abs(nodos(:,2) - L/2) <1e-6; %y = L/2

mesh.cond_borde(borde_CD,[1,3,5]);
mesh.cond_borde(borde_AB,sym_xz);
mesh.cond_borde(borde_AC,sym_yz);

%% Cargas

mesh.armar_R(q);

%% Calculos

mesh.calc_U;

%% Graficar
mesh.dibujar;

%% Valores del benchmark
wB = mesh.U(mesh.dofs(mesh.counts.nnod,3)) %el ultimo nodo en z