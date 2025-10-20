clc
clear

L = 10;

q = 10000;%N/m
elems = 5;
x = linspace(0,L,elems+1)';
y = zeros(elems+1,1);
mesh.nodos = [x y];

x = 1:elems;
y = 2:elems+1;
mesh.elems.con = [x;y]';

mesh.elems.E = 210e9;%Gpa
mesh.elems.I = 10e-4;

nnod = size(mesh.nodos,1);
nelem = size(mesh.elems.con,1);
ndof = nnod*2;

mesh.dofs = reshape(1:ndof,2,[])';

mesh.free = true(nnod,2);
mesh.free(1,:) = false;
mesh.free(end,1) = false;
mesh.free = reshape(mesh.free',1,[]);

R = zeros(ndof,1);
L = L/nelem;
for i = nelem
    R(2*i-1) = -q*L/2;
    R(2*i) = -q*L^2/12;
    R(2*i+1) = -q*L/2;
    R(2*i+2) = q*L^2/12;
end

K = zeros(ndof);
for i = 1:nelem
    nodoid = mesh.elems.con(i,:);

    coord = mesh.nodos(nodoid,:);
    V = coord(2,:)-coord(1,:);
    L = norm(V);

    E = mesh.elems.E;
    I = mesh.elems.I;
    
    Y1 = 12*E*I/L^3;
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kloc = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];
    
    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0 0;
         0 0 c s];
    
    Kglob = Kloc;

    dir = mesh.dofs(nodoid,:)';

    K(dir,dir) = K(dir,dir) + Kglob;
end

Kr = K(mesh.free,mesh.free);
Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

F = K*U;

syms x Li
X = [1 x x^2 x^3];
A = [subs(X, x ,0)
     subs(diff(X,x), x, 0)
     subs(X, x, Li)
     subs(diff(X,x), x, Li)];

N = X/A;

divisiones = 6;
def = [];
Largo = [];
Coord = mesh.nodos;
Elem = mesh.elems.con;
for i = 1:nelem
    V = (Coord(Elem(i, 2),:) - Coord(Elem(i, 1),:));
    Le = norm(V);

    dir = [Elem(i,1)*2-1 Elem(i,1)*2 Elem(i,2)*2-1 Elem(i,2)*2];
    
    Uloc = U(dir);

    Vy = subs(N, Li, Le)*Uloc;   

    Inicio = Coord(Elem(i,1), 1);

    for j = 0:Le/divisiones:Le
        def = [def subs(Vy, x, j)];
        Largo = [Largo Inicio+j];
    end
end

plot(Largo, def);
