clc
clear

P = 1000;
L = 5;
M = 1;

E = 200e9;
v = 0.3;
t = 0.05;

C = t*E/(1-v^2)*[1 v 0;
                v 1 0;
                0 0 (1-v)/2];


nelem = 10000;

%saltos entre filas de nodos
%33-34...52
%21-22-..32
%1-2-...20
salto = nelem*2+1;
salto2 = salto+nelem+1;
padding = 0:2;

%distancias entre nodos
x1 = linspace(0,L,salto);%arriba y abajo del elemento
x2 = linspace(0,L,nelem+1);%los 2 del medio
y = linspace(0, M, 3);%distancia en y

%creo los puntos que dividen el espacio y los elementos entre los nodos
for i = 1:salto
    mesh.nodes(i,:) = [x1(i) y(1)];
    mesh.nodes(salto2+i:salto2+i,:) = [x1(i) y(3)];
end
for i = 1:nelem+1
    mesh.nodes(salto+i,:) = [x2(i) y(2)];
end

for i = 1:nelem
    primero = 2*i-1;
    mesh.elem(i,:) = [primero+padding primero+salto+1-i primero+salto+2-i primero+salto2+padding];
    mesh.elem(i,2:end) = mesh.elem(i,[3 8 6 2 5 7 4]);%ponerlo como dice el crear K
end

nnod = size(mesh.nodes,1);
ndof = nnod*2;

mesh.dofs = reshape(1:ndof, 2, [])';

mesh.Libres = true(nnod, 2);
empotrados = [1 salto+1 salto2+1];%empotro todos los de la izquierda
mesh.Libres(empotrados,:) = false;

mesh.Libres = reshape(mesh.Libres', [], 1);

R = zeros(ndof,1);
R(2*salto-1) = P;
R(2*(salto+salto2)-1) = -P;
K = spalloc(ndof,ndof,32*ndof);

for i = 1:nelem
    dir_nod = mesh.elem(i,:);

    coord_nodos = mesh.nodes(dir_nod,:);
    
    Kloc = crearK_Q8(coord_nodos, C);

    dir = mesh.dofs(dir_nod,:)';

    K(dir,dir) = K(dir,dir) + Kloc;
end

Rr = R(mesh.Libres);

Kr = K(mesh.Libres, mesh.Libre

U = zeros(ndof,1);
U(mesh.Libres) = Kr\Rr;

divisiones = 20;
mult = 10000; %multiplica los desplazamientos

plotQ8(mesh.nodes,mesh.elem, mesh.dofs, U, divisiones, mult)