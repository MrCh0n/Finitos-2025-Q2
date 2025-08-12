clc
clear

nelem = 10;

Largo = 40; %m
Peso = 110000*9.81; %N

b = 0.250;
h = 1;
t = 4*25.4/1000;

I = (b*h^3 - (b-2*t)*(h-2*t)^3)/12;

q = Peso/Largo;
L = Largo/nelem;

E = 200e9;

mesh.coord = (0:nelem)*L;

nnod = size(mesh.coord,2);
ndof = nnod*2;

for i = 1:nelem
    mesh.elem(i,:) = [i i+1];
end

K = zeros(ndof);

mesh.dofs = 1:ndof;
mesh.dofs = reshape(mesh.dofs,2,[])';

Libres = true(nnod,2);

Libres([1 nnod], 1) = false;%simplemetne apoyado con 1, poner : para doble empotrado

Libres = reshape(Libres', 1, []);

R = zeros(ndof,1);

for i = 1:nelem
    %ya tengo el L de antes
    Y1 = 12*E*I/L^3;
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kel = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];

    dirnod1 = mesh.dofs(mesh.elem(i,1) , :);
    dirnod2 = mesh.dofs(mesh.elem(i,2) , :);
    
    P = q*L/2;
    M = q*L^2/12;
    
    dir = [dirnod1 dirnod2];

    R(dir) = R(dir) + [-P; -M; -P; M];

    K(dir, dir) = K(dir, dir) + Kel;
end

Rr = R(Libres);

Kr = K(Libres, Libres);

U = zeros(ndof,1);
U(Libres) = Kr\Rr;

%divisiones por elemento
divisiones = 5;

stress = zeros(divisiones*nelem,1);
for i = 1:nelem
    x = linspace(0, L, divisiones);
    B = [-6/L^2+12*x/L^3; -4/L+6*x/L^2; 6/L^2-12*x/L^3; -2/L+6*x/L^2]';

    dir = [mesh.dofs(mesh.elem(i,1) , :) mesh.dofs(mesh.elem(i,2) , :)];

    stress((i-1)*divisiones+1:i*divisiones) = h/2*E*B*U(dir);
end

Le = L;
%TODO hacer el B
syms x L
X = [1 x x^2 x^3];
A = [subs(X, x ,0)
     subs(diff(X,x), x, 0)
     subs(X, x, L)
     subs(diff(X,x), x, L)];

N = X/A;

%TODO las tensiones
def = [];
Largo = [];

for i = 1:nelem
    dir = [mesh.dofs(mesh.elem(i,1) , :) mesh.dofs(mesh.elem(i,2) , :)];

    Uloc = U(dir);

    Vy = subs(N, L, Le)*Uloc;   

    Inicio = mesh.coord(mesh.elem(i,1));

    for i = 0:Le/divisiones:Le
        def = [def subs(Vy, x, i)];
        Largo = [Largo Inicio+i];
    end
end

plot(Largo, def);
