% Main del trabajo practico 1

clc
clear all
close all

addpath(genpath('C:\Users\franc\OneDrive - ITBA\ITBA\ITBA\9no Cuatrimestre\FEM\Codigos Matlab\Chon\Finitos-2025-Q2\Libreria_elementos'))

%% Datos del problema
E = 4.32e8; %psi
v = 0;

R = 25; %in
L = 25; %in
t = 0.25 ; %in
tita = 40; %grados

q = -90; %psi

%% Control
type = 1; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 4; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        T = t;
        dofselem = 6;
        crearK = @crearK_shellMQ4;
        sym =
    case 2
        T = t*ones(1,4);
        dofselem = 5;
        crearK = @crearK_shell_degeneradoQ4;
        sym = 
end
%% Malla

[nodos,elems] = mallador_ej1(L/2,R,tita,div,div);

figure(2)
draw_Mesh(elems,nodos, 'NodeLabel',true,'Type','Q4','Color','b')
hold off


nnod = size(nodos,1); %cant nodos
ndof = dofselem*nnod; 
nelem = size(elems,1); %cant elementos

dofs = reshape(1:ndof,dofselem,[])';

%% Calcular direccion z

%% Matriz de rigidez

K = zeros(ndof);


for i=1:nelem    
    nodoid = elems(i,:);

    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    
    Kel = crearK(nodos(nodoid,:),E,v,T);

    K(dir,dir)=K(dir,dir) + Kel;
end

%% Condiciones de borde
free = true(ndof,1);

borde_CD = abs(nodos(:,1)) <1e-6; %x = 0
borde_AC = abs(nodos(:,2)) <1e-6; %y = 0
borde_AB = abs(nodos(:,2) - L/2) <1e-6; %y = L/2


free(dofs(borde_CD,:)) = false; %empotrado
free(dofs(borde_AC,[1 5 6])) = false; %sym (mov en x, giro en y z)
free(dofs(borde_AB,[2 4 6])) = false; %sym (mov y, giro en x y z)

%% Cargas

% Definición de las cargas aplicadas
R = zeros(ndof, 1);
for i=1:nelem
    nodoid = elems(i,:);
    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    Ae = area(nodos(nodoid,:));
    

    R(dir(3:dofselem:4*dofselem)) = R(dir) + Ae*q/4; % Carga en z
end


%% Calculos

Kr = K(free,free);

Rr = R(free);

%Calculo movimientos
U = zeros(ndof,1);
U(free) = Kr\Rr;

%% Graficar

%% Funciones


function [Ae] = area(nodos)
    v1 = (nodos(2,:)-nodos(1,:))/norm(nodos(2,:)-nodos(1,:));
    
    v3 = cross(v1,(nodos(4,:)-nodos(1,:)));
    
    v3 = v3 / norm(v3); % Normalize the cross product vector
    
    v2 = cross(v3,v1);
    
    Q = [v1;v2;v3]; %v1 es vector fila
    nodosp = (Q*nodos')'; %nodos locales
    nodos = nodosp(:,1:2);

    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];
    [w, puntos, n] = gauss([2,2]);
    Ae = 0; % area del elemento
    for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        mult = abs(det(J))*w(i);
        Ae = Ae + mult;
    end
end