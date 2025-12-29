clc
clear

k = 1;
t = 1;

R = 6.37e6;

n = 100;%cantidad de nodos en x
m = 100;%cantidad de nodos en y

L = 10;
W = 10;

G = 6.674e-11;
Mt = 5.97e24;

M = Mt*G*4*pi;

divx = n-1;
divy = m-1;
bordes = [0 0; L 0; L W; 0 W];

[nodos, elems] = mallador_cuadrado_Q4(bordes, divx, divy);

nnod = size(nodos,1);
ndof = nnod;
nelem = size(elems,1);

izquierda = 1:divy+1;
arriba = divy+1:divy+1:nnod;
abajo = 1:divy+1:nnod;
derecha = nnod-divy:nnod;

K = spalloc(ndof,ndof, 32*ndof);

dofs = 1:ndof;

Q = zeros(ndof,1);

masa = elems(ceil(divy*divx/2),:);

elem1 = ceil(divy*divx/2)+floor(divy/2.2);
elem2 = ceil(divy*divx/2)-floor(divy/2.2);
elems(elem2,:)
elems(elem1,:)
coord = [elems(elem2,1:2) elems(elem1,1:2)];

for i = -5:5
    el1 = elem1+2*i*divy;
    el2 = elem2+2*i*divy;

    coord = [elems(el2,1:2) elems(el1,1:2)];

    elems(el2-divy,2) = coord(3);
    elems(el2+divy,1) = coord(4);
    elems(el2,1:2) = coord(3:4);

    elems(el1-divy,2) = coord(1);
    elems(el1+divy,1) = coord(2);
    elems(el1,1:2) = coord(1:2);
end
%Q(masa) = M/4;

free = true(nnod, 1);

fijos = [izquierda derecha abajo arriba];

free(fijos) = false;

for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = crearK_termico_Q4(coord,k);

    dir = dofs(nodoid);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
        
    K(dir,dir) = K(dir,dir) + Kel;
end

T = zeros(ndof,1);

T(abajo) = -G*Mt/R;
T(arriba) = -G*Mt/(R+W);
up = linspace(R,R+W,m);
T(derecha) = -G*Mt./up;
T(izquierda) = -G*Mt./up;
T(fijos) = T(fijos) +G*Mt/R;

Q = Q + K*T;

Kr = K(free,free);

Qr = Q(free);

T(free) = -Kr\Qr;


X = nodos(:,1);
Y = nodos(:,2);
Tim = reshape(T,m,n);

imagesc(X,Y,Tim)
colorbar

figure(2)

U = zeros(nelem,1);
V = zeros(nelem,1);
cont = zeros(nelem,1);

for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    dir = dofs(nodoid);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

    Tel = T(dir);
    
    [Uel, Vel] = flujo_Q4(coord,Tel,k);
    U(i) = U(i) + Uel;
    V(i) = V(i) + Vel;
    cont(i) = cont(i) + 1;
end

U = U./cont;
V = V./cont;

% x = linspace(0,L,n);
% y = linspace(0,W,m);
% [X,Y] = meshgrid(x,y);
% contour(X,Y,Tim,20)
% colorbar

hold on
x = linspace(0,L,divx);
y = linspace(0,W,divy);
[X,Y] = meshgrid(x,y);
X = reshape(X,1,[])';
Y = reshape(Y,1,[])';

quiver(X,Y,U,V);
axis equal
hold off

figure(3)

gravedad = (U.^2+V.^2).^0.5;

gravedad = reshape(gravedad,divx,divy);

imagesc(gravedad)

colorbar