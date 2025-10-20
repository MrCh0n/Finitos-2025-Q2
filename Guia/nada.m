clc
clear

caso = 1;
dt = 25;

P = 1000;

alpha_acero = 12e-6;
alpha_bronce = 17e-6;

mesh.mat(1).E = 210*1e9;
mesh.mat(2).E = 100*1e9;

De = 100e-3;
Di = 0;

I = pi*(De^4-Di^4)/64;
A = pi*(De^2-Di^2)/4;
y_max = De/2;

mesh.nodos = [0 0;
              1 0
              2 0];

mesh.elems.con = [1 2;
                  2 3];

nnod = size(mesh.nodos,1);
ndof = nnod*3;
nelem = size(mesh.elems.con,1);

mat = ones(nelem,1);

mesh.dofs = reshape(1:ndof,3,[])';

mesh.free = true(nnod,3);
empotrados = 1;
mesh.free(empotrados,:) = false;

R = zeros(ndof,1);

R(end-2) = P;

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

for i = 1:nelem 
    [sigma, Ft] = cargas_termicas(A, mesh.mat(i).E, alpha_acero, dt);%barra acero 
    puesto = 3*i-2;
    R([puesto puesto+3]) =  R([puesto puesto+3]) + Ft;
    s_t(i) = sigma;
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