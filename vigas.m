function [] = vigas(sistema)

 campos = fieldnames(sistema);

Coord = sistema.(campos{1});
Elem = sistema.(campos{2});
BC = sistema.(campos{3});
Fnodos = sistema.(campos{4});


nnod = size(Coord,1); %cant de nodos
ndof = 2*nnod; %degrees of freedom
nelem = size(Elem,1); %cantidad de elementos
K = zeros(ndof);

%armar la matriz K
for i = 1:nelem
    V = (Coord(Elem(i, 2),:) - Coord(Elem(i, 1),:));
    L = norm(V);

    E = Elem(i,3);
    I = Elem(i,4);
    
    Y1 = 12*E*I/L^3;
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kel = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];

    dir = [Elem(i,1)*2-1 Elem(i,1)*2 Elem(i,2)*2-1 Elem(i,2)*2];
    
    %TODO poner la rotacion
    K(dir,dir) = K(dir,dir) + Kel;
end

Libres = 1:ndof;
Libres(BC) = [];

Kreducido = K(Libres, Libres);

R = zeros(ndof,1);
for i = 1:size(Fnodos, 1)
    R(Fnodos(i,1)) = Fnodos(i,2);
end

Rreducido = R(Libres);

U = zeros(ndof, 1);
U(Libres) = Kreducido\Rreducido

%TODO hacer el B
syms x L
X = [1 x x^2 x^3];
A = [subs(X, x ,0)
     subs(diff(X,x), x, 0)
     subs(X, x, L)
     subs(diff(X,x), x, L)];

N = X/A;

%TODO las tensiones

divisiones = 6;
def = [];
Largo = [];

for i = 1:nelem
    V = (Coord(Elem(i, 2),:) - Coord(Elem(i, 1),:));
    Le = norm(V);

    dir = [Elem(i,1)*2-1 Elem(i,1)*2 Elem(i,2)*2-1 Elem(i,2)*2];

    Uloc = U(dir);

    Vy = subs(N, L, Le)*Uloc;   

    Inicio = Coord(Elem(i,1), 1);

    for j = 0:Le/divisiones:Le
        def = [def subs(Vy, x, j)];
        Largo = [Largo Inicio+j];
    end
end

plot(Largo, def);
