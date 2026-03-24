clc
%clear

x1 = [-1; 1; 1; -1; -1; 1; 1; -1];
y1 = [-1; -1; 1; 1; -1; -1; 1; 1];
z1 = [-1; -1; -1; -1; 1; 1; 1; 1];

nodos = [x1 y1 z1];

elem = [1 2 3 4 5 6 7 8];

nnod = size(nodos,1);
ndof = nnod*3;
nelem = size(elem,1);

dofs = 1:ndof;
dofs = reshape(dofs,3,[])';

empotrados = [1 4 5 8];
free = true(nnod, 3);

free(empotrados,:) = false;

free = reshape(free',1,[]);

E = 200e9;
v = 0.3;

A = (1-v)*E/((1+v)*(1-2*v));
B = v*E/((1+v)*(1-2*v));
G = E/(2*1+v);

C = [A B B 0 0 0;
    B A B 0 0 0;
    B B A 0 0 0;
    0 0 0 G 0 0;
    0 0 0 0 G 0;
    0 0 0 0 0 G];

K = zeros(ndof);
R = zeros(ndof,1);

carga_x = [2 3 6 7]*3-2;
carga_x2 = [6 7]*3-2;


R(carga_x) = 1;
%R(carga_x2) = -1;

for i = 1:nelem
    nodoid = elem(i,:);

    coord = nodos(nodoid,:);

    eledofs = reshape(dofs(nodoid,:)',1,[]);

    Kel = crearK_H8(coord,C);

    K(eledofs,eledofs) = K(eledofs,eledofs) + Kel;
end

Kr = K(free,free);
Rr = R(free);

U = zeros(ndof,1);
U(free) = Kr\Rr;

mov = reshape(U,3,[])';


