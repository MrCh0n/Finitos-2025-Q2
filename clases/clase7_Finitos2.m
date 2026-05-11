clc
clear
close all

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos
L = 1;%m
W = 0.1;%m

E = 10e9;%Gpa
v = 0.3;
t = 0.1;
alpha = 1;
Pp = 0.5e6;%Mpa

F = 0;%1e6*L*t;% Carga total [N]

%Stress
C = E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
Czz = zeros(1,3);

bordes = [0 0;
          L 0;
          L W;
          0 W];

divx = 40;
divy = 2;
%% Mesh
[nodos, elems, bordes] = mallador_cuadrado_Q4(bordes, divx, divy);

nnod = size(nodos,1);
ndof = nnod*2;
nelem = size(elems,1);

dofs = reshape(1:ndof,2,[])';

free = true(nnod, 2);
%empotrado a izquierda
empotrado = bordes.lado_41;
free(empotrado,:) = false;

free = reshape(free', 1, []);
%% K y Cg
dofselem = 8;

nnz = nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros

I = zeros(nnz,1);
J = zeros(nnz,1);
V = zeros(nnz,1);
cont = 1;

Cg = zeros(ndof, nnod);
for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = t*crearK_Q4(coord, C);

    Cg_el = crearCg_Q4(coord, t, alpha);

    dir = reshape(dofs(nodoid,:)',1,[]);

    Cg(dir,nodoid) = Cg(dir,nodoid) + Cg_el;
    for a = 1:dofselem
        for b = 1:dofselem
            I(cont) = dir(a);
            J(cont) = dir(b);
            V(cont) = Kel(a,b);
            cont = cont+1;
        end%b
    end%a
end

K = sparse(I,J,V);
%% P y R
P = ones(nnod,1)*Pp;

R = zeros(ndof,1);

q = F/divx;% la fuerza en los nodos centrales
arriba = bordes.lado_34;
arr_y = dofs(arriba,2);

R(arr_y) = R(arr_y) - q;
R(arr_y([1 end])) = R(arr_y([1 end]))/2;

R = R - Cg*P;
%% Calculo
U = zeros(ndof,1);

Kr = K(free,free);
Rr = R(free);

U(free) = Kr\Rr;

w = F/L;
I = W^3*t/12;
delta = -w*L^4/(8*E*I);
error = (min(U)-delta)/abs(delta)*100

%Stress
[bruto, n_bruto, n_suave] = stress_2D(nodos, elems, dofs, U, C, Czz, @stress_Q4, @global_Q4, @elem_a_nodos_Q4);
%% Plot
escala = 1;

x = nodos(:,1);
y = nodos(:,2);

x_deformada = x + escala*U(1:2:ndof);
y_deformada = y + escala*U(2:2:ndof);
nodos_deformada = [x_deformada y_deformada];

figure()
draw_Mesh(elems, nodos,'Type','Q4','Color','b')
hold on
draw_Mesh(elems,nodos_deformada,'Type','Q4','Color','k')
hold off

%draw_stress(nodos,elems,bruto)