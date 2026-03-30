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
type = 1; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 16; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        mesh = Mindlin(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5 6];
        sym_xz = [2 4 6];
    case 2
        T = t*ones(1,4);
        dofselem = 5;
        crearK = @crearK_shell_degeneradoQ4;
        sym_yz = [1 5];
        sym_xz = [2 4];
end
%% Malla
%ya esta mallado
nodos = mesh.nodes;
elems = mesh.elems;
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

function [Ae] = Volumen_degenerado(nodos,T)
    % T es 4x1 de los espesores en el elemento
    % nodos es 4x3 de las coordenadas
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];
    [w, puntos, n] = gauss([2,2]);
    Ae = 0; % area del elemento
    for i = 1:n
        N = [1 puntos(i,1) puntos(i,2) puntos(i,1)*puntos(i,2)]/A;
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];     
        J = D * nodos;         % 2x3 jacobiano de xy->xyz
        
        X = J(1,:);          
        Y = J(2,:);        

        J = [norm(X) 0;
             0        norm(Y)];

        espesor = N*T;
        mult = abs(det(J))*w(i);
        Ae = Ae + mult*espesor;
    end
end