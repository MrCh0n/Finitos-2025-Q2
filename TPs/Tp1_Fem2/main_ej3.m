% Main del trabajo practico 1

clc
clear all
close all

%% Datos del problema
E = 3.1e6; %psi
v = 0.3;

R = 300; %in
L = 600; %in
t = 3 ; %in
tita = 90; %grados

q = -90; %psi
P = 1;%lbf

%% Control
type = 2; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 50; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        dofselem = 6;
        mesh = Mindlin(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5 6];
        sym_xz = [2 4 6];
        sym_xy = [3 4 5];
    case 2
        dofselem = 5;
        mesh = Degenerado(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5];
        sym_xz = [2 4];
        sym_xy = [3 4 5];
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
borde_BD = abs(nodos(:,3)) < 1e-6;

mesh.cond_borde(borde_CD,[1,3,5]);
mesh.cond_borde(borde_AB,sym_xz);
mesh.cond_borde(borde_AC,sym_yz);
mesh.cond_borde(borde_BD,sym_xy);

%% Cargas

%mesh.armar_R(q);
R = zeros(mesh.counts.ndof, 1);
cargas = mesh.dofs(div+1,3);
R(cargas) = -P;

mesh.cargar_R(R);


%% Calculos

mesh.calc_U;

%% Graficar
mesh.dibujar;

%% Valores del benchmark
wB = mesh.U(cargas) %el ultimo nodo en z