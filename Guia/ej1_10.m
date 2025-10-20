clc
clear

caso = 1;
dt = 25;

alpha_acero = 12e-6;
alpha_bronce = 17e-6;

mesh.mat(1).E = 210*1e9;
mesh.mat(2).E = 100*1e9;
mat = [1 1 1 2];

De = 100e-3;
Di = 0;

I = pi*(De^4-Di^4)/64;
A = pi*(De^2-Di^2)/4;
y_max = De/2;

mesh.nodos = [0 0;
              0 250;
              200 250;
              280 250;%Barra bronce
              280 0]*1e-3;%barra broce

mesh.elems.con = [1 2;
                  2 3;
                  3 4;
                  4 5];

nnod = size(mesh.nodos,1);
ndof = nnod*3;
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,3,[])';

mesh.free = true(nnod,3);
empotrados = [1 5];
mesh.free(empotrados,:) = false;
mesh.free(3,[1 2]) = false;%fijo

R = zeros(ndof,1);

K = zeros(ndof);

for i =1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);
    
    E = mesh.mat(mat(i)).E;

    Kel = crearK_viga_barra(nodos,E,A,I);

    dofs1 = mesh.dofs(mesh.elems.con(i,1),:);
    dofs2 = mesh.dofs(mesh.elems.con(i,2),:);
    eleDofs = [dofs1 dofs2];%donde tiene que ir

    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Kel;
end

U = zeros(ndof,1);
s_t = zeros(nelem,1);
if caso == 1 || caso == 4
    [sigma, Ft] = cargas_termicas(A, mesh.mat(2).E, alpha_acero, dt);%barra acero 
    R([2 5]) = R([2 5]) + Ft;
    s_t(1) = sigma;
end

if caso == 2 || caso == 4
    [sigma, Ft] = cargas_termicas(A, mesh.mat(2).E, alpha_acero, dt);%barra acero
    R([4 7]) = R([4 7]) + Ft;
    s_t(2) = sigma;

    R([7 10]) = R([7 10]) + Ft;
    s_t(3) = sigma;
end

if caso == 3 || caso == 4
    [sigma, Ft] = cargas_termicas(A, mesh.mat(2).E, alpha_bronce, dt);%barra bronce   
    R([11 14]) = R([11 14]) + -Ft;
    s_t(4) = sigma;
end

mesh.free = reshape(mesh.free',[],1);

Kr = K(mesh.free,mesh.free);
Rr = R(mesh.free);


U(mesh.free) = Kr\Rr;

div = 10;

stress = zeros(nelem,5);
minses = zeros(nelem,div);
for i = 1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);

    E = mesh.mat(mat(i)).E;

    dofs1 = mesh.dofs(mesh.elems.con(i,1),:);
    dofs2 = mesh.dofs(mesh.elems.con(i,2),:);
    eleDofs = [dofs1 dofs2];%donde tiene que ir
    
    Uel = U(eleDofs);

    [stress(i,:), minses(i,:)] = stress_viga_barra(nodos,Uel,E,y_max,div,s_t(i));
end