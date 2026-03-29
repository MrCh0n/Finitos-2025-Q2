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

q = 90; %psi

%% Control
type = 1; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 4; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        T = t;
        dofselem = 6;
        crearK = @crearK_shellMQ4;
    case 2
        T = t*ones(1,4);
        dofselem = 5;
        crearK = @crearK_shell_degeneradoQ4;
end
%% Malla
bordes = [0,0;
          L,0;
          L,L/2;
          0,L/2];
[nodos,elems] = mallador_ej1(bordes,div,div,R,tita);

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
%borde_BD = abs(nodos(:,1) - a) <1e-6;
%borde_AB = abs(nodos(:,2) - a) <1e-6;

% empotrado
%dofs2fix = [dofs(borde_1,:);dofs(borde_2,:);dofs(borde_3,:);dofs(borde_4,:)];

free(dofs2fix) = false;