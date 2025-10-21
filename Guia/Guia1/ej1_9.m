clc
clear
close all

auxiliar = [0 1 0];%direccion

D = 40e-3;%m diametro
%Tubo hueco y acero
E = 210e9;
G = 70e9;
rho = 7850;
Iz = pi*D^4/64;
Iy = Iz;
A = pi*D^2/4;
Ip = Iz*2;
y_max = D/2;%m
z_max = D/2;%m
r = D/2;%tiene sentido para un IPN?

nodos = [0 0 0;%A
         5 0 0;%E
         10 0 0;%B
         20 0 0;%C
         26 0 0;%F
         32 0 0]*25.4e-3;%D

elems = [1 2;
        2 3;
        3 4;
        4 5;
        5 6];

nelem = size(elems,1);

div = ones(nelem,1);

[mesh.nodos, mesh.elems.con, mesh.elems.transformada, mesh.elems.inicio] = mesh_1D(elems, nodos,div);

nnod = size(mesh.nodos,1);
ndof = 6*nnod;%viga 3D
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,6,[])';

mesh.free = true(nnod,6);
rotacion_en_x = mesh.elems.transformada([1 3 4 6]);%A B C D
punto_A = mesh.elems.transformada(1);
punto_E = mesh.elems.transformada(2);

mesh.free(rotacion_en_x,[2 3 5 6]) = false;
mesh.free(punto_A,1) = false;%agarra en x
mesh.free(punto_E,4) = false;%no gira
mesh.free = reshape(mesh.free',[],1);

cargaF = [0 -3500 0 880 0 0];%N
cargaE = [2000 -1600 0 0 0 0];%N

R = zeros(ndof,1);
punto_E = mesh.dofs(mesh.elems.transformada(2),:);
punto_F = mesh.dofs(mesh.elems.transformada(5),:);

R(punto_E) = cargaE;
R(punto_F) = cargaF;

K = zeros(ndof);
for i = 1:size(div,1)
    nodoid = mesh.elems.con(mesh.elems.inicio(i), :);

    nodos = mesh.nodos(nodoid,:);

    Kglobal = crearK_viga_3D(nodos,E,G,A,Iy,Iz,Ip, auxiliar);

        %meterlo en la matriz enorme
    cantidad = 0:div(i)-1;
    for j = (mesh.elems.inicio(i) + cantidad)
        dofs1 = mesh.dofs(mesh.elems.con(j,1),:);
        dofs2 = mesh.dofs(mesh.elems.con(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        
        K(eleDofs,eleDofs) = K(eleDofs, eleDofs) + Kglobal;
        %KelGlobal; 
        %spy(K);
    end 
end

Kr = K(mesh.free,mesh.free);
Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

F = K*U - R;

girosx =  U(4:6:end);
giro_A_D = girosx(mesh.elems.transformada(6)) - girosx(mesh.elems.transformada(1));%radianes
giro_A_D = giro_A_D*180/pi 

%% Stress y graficado
divisiones = 5;%divisiones por subelemto

stress = zeros(nelem,5);
for i = 1:size(div,1)
    nodoid = mesh.elems.con(mesh.elems.inicio(i), :);

    nodos = mesh.nodos(nodoid,:);

        %meterlo en la matriz enorme
    cantidad = 0:div(i)-1;
    for j = (mesh.elems.inicio(i) + cantidad)
        dofs1 = mesh.dofs(mesh.elems.con(j,1),:);
        dofs2 = mesh.dofs(mesh.elems.con(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        
        Uel = U(eleDofs);

        stress(j,:) = stress_viga_3d(nodos, Uel, E, G, y_max, z_max, r, auxiliar, divisiones);
        %KelGlobal; 
        %spy(K);
    end 
end

elementos = [mesh.elems.inicio; nelem];
mult = 5;
for i = 1:size(div,1)
    nodoid = mesh.elems.con(elementos(i):elementos(i+1), 1);
    nodos = mesh.nodos(nodoid,:);
    
    eledofs = mesh.dofs(nodoid,:);
    eledofs = reshape(eledofs', [], 1);

    Uel = U(eledofs);
    draw_viga_3D(nodos, Uel,mult)
end
