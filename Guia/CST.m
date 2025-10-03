clc
clear

a = 5;

n = 3*a+1;%tiene que ser multiplo de 3 + 1
m = 2*a+1;%para que quede cuadrado

P = 2000; %N 

v = 0.3; %poisson acero
E = 200e9; %Pa
t = 2/1e3; %m

Lx = 6; %m
Ly = 4; %m

x = linspace(0, Lx, n);
y = linspace(0, Ly, m);

for i = 1:m
    indice = (1:n);
    indice = indice + (i-1)*n;
    mesh.coord(indice,1) = x;
    mesh.coord(indice,2) = y(i);
end

%necesito hacer 1 menos de la cantidad de nodos
%y como los quiero usar para hacer indice le tengo que sacar 1 cada n
for i = 0:m-2
    inicio = i*n;
    for j = 1:n-1
        dir = inicio+j-i;
        mesh.elem(2*(dir)-1,:) = inicio + [j j+1 j+n];
        mesh.elem(2*(dir),:) = inicio +[j+1 j+n j+n+1];% crea de a 2 triangulos
    end
end

nnod = size(mesh.coord,1);
ndof = 2*nnod;
nelem = size(mesh.elem,1); 

mesh.dofs = reshape(1:ndof, 2, [])';

mesh.free = true(nnod,2);

emp = (n*(m-1)+1);
mesh.free(emp:-n:emp-n*(a-1), 2) = false;%empotramiento

for i = 1:n:n*m
    mesh.free(i, :) = false;%le saco X a todos los que estan en la izquierda
end

mesh.free = reshape(mesh.free', [], 1);

R = zeros(ndof,1);
carga_dir = (n-1)*[2 4 6]/Lx + 1;%pone la carga a 2 4 6 m de la zona inferior
R(2*carga_dir) = -P;

C = E*[1 v 0;
       v 1 0;
       0 0 (1-v)/2]/(1-v^2); %constitutiva

K = zeros(ndof);

for i = 1:nelem
    dir_nod = mesh.elem(i, :);

    nodos = mesh.coord(dir_nod,:);
    
    Kloc = t*crearK_CST(nodos, C);

    dir = mesh.dofs(dir_nod,:)';
    
    K(dir, dir) = K(dir, dir) + Kloc;
end

Kr = K(mesh.free, mesh.free);
Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

mult = 1000;

plot_CST(mesh.coord, mesh.elem , mesh.dofs, U, mult, C)