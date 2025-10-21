clc
clear
close all

carga = [0 -10 0 0 0 0]*1e3;%Kn

R = zeros(ndof,1);
pos = mesh.dofs(mesh.elems.transformada(5),:);
R(pos) = carga;

auxiliar = [0 1 0];%direccion

%IPN200 y acero
E = 210e9;
G = 70e9;
Iz = 2140e-8;
Iy = 117e-8;
A = 33.4e-4;
Ip = 1e-8;
y_max = 200e-3/2;%m
z_max = 90e-3/2;%m
r = (y_max+z_max)/2;%tiene sentido para un IPN?

nodos = [4 0 0;
         4 0 -6;
         0 0 -6;
         0 0 0;
         4 0 -3];

elems = [4 1;
        1 5;
        5 2;
        2 3;
        3 4];

nelem = size(elems,1);

div = ones(nelem,1);

[mesh.nodos, mesh.elems.con, mesh.elems.transformada, mesh.elems.inicio] = mesh_1D(elems, nodos,div);

nnod = size(mesh.nodos,1);
ndof = 6*nnod;%viga 3D
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,6,[])';

mesh.free = true(nnod,6);
empotrados = mesh.elems.transformada(3:4);

mesh.free(empotrados,:) = false;
mesh.free = reshape(mesh.free',[],1);

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
