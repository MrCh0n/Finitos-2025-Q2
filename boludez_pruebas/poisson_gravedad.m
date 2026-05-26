clc
clear

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos
k = 1;%vale 1 solo porque lo necesito para el crearK y flujo

n = 300;%cantidad de nodos en x
m = 300;%cantidad de nodos en y

L = 10;%ancho de la ventana
W = 10;%alto de la ventana

x0 = L/4;%centro de la masa
y0 = 3*W/4;%centro de la masa
p = 0;%masa en el centro de la masa

G = 6.674e-11;% Constante gravitacional
Mt = 5.97e24;% Masa tierra
R = 6.37e6;% Radio tierra

M = Mt*G*4*pi;% poisson de la tierra

%% Masheado
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

dofs = 1:ndof;

Q = zeros(ndof,1);

masa = elems(ceil(divy*divx/2),:);

elem1 = ceil(divy*divx/2)+floor(divy/2.2);
elem2 = ceil(divy*divx/2)-floor(divy/2.2);
elems(elem2,:);
elems(elem1,:);
coord = [elems(elem2,3:4) elems(elem1,1:2)];
cant_elems = ceil(n/20);
%Q(masa) = M/4;

%draw_Mesh(elems,nodos,"Type","Q4","NodeLabel",true)

free = true(nnod, 1);

fijos = [izquierda derecha abajo arriba];

free(fijos) = false;

elem_movido = elems;
 for i = -cant_elems:cant_elems
     el1 = elem1+2*i*divy;
     el2 = elem2+2*i*divy;

     coord = [elems(el2,1:2) elems(el1,1:2)];

    elem_movido(el2-divy,2) = coord(3);
    elem_movido(el2+divy,1) = coord(4);
    elem_movido(el2,1:2) = coord(3:4);

    elem_movido(el1-divy,2) = coord(1);
    elem_movido(el1+divy,1) = coord(2);
    elem_movido(el1,1:2) = coord(1:2);
 end

dofselem = 4;%cantidad de dofs totales por elemento
nnz = nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros

I = zeros(nnz,1);
J = zeros(nnz,1);
V = zeros(nnz,1);

cont = 1;

for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = crearK_termico_Q4(coord,k);
    
    nodoid = elem_movido(i,:);
    dir = dofs(nodoid);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

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

T = zeros(ndof,1);

T(abajo) = -G*Mt/R;
T(arriba) = -G*Mt/(R+W);
up = linspace(R,R+W,m);
T(derecha) = -G*Mt./up;
T(izquierda) = -G*Mt./up;
T(fijos) = T(fijos) +G*Mt/R;

Q = Q + K*T;
for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    den = ((coord(:,1)-x0).^2+(coord(:,2)-y0).^2+1);
    rho = p./den;

    dir = dofs(nodoid);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
        
    Q(dir) = Q(dir) + rho;
end

Kr = K(free,free);

Qr = sparse(Q(free));

T(free) = -Kr\Qr;

%% Ploteo
X = nodos(:,1);
Y = nodos(:,2);
Tim = reshape(T,m,n);
figure(4)
imagesc(X,Y,Tim)
colorbar

figure(2)

U = zeros(nelem,1);
V = zeros(nelem,1);
cont = zeros(nelem,1);
xc = zeros(nelem,1);
yc = zeros(nelem,1);

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
    xc(i) = mean(coord(:,1));
    yc(i) = mean(coord(:,2));
end

U = U./cont;
V = V./cont;

% x = linspace(0,L,n);
% y = linspace(0,W,m);
% [X,Y] = meshgrid(x,y);
% contour(X,Y,Tim,20)
% colorbar

hold on
% x = linspace(0,L,divx);
% y = linspace(0,W,divy);
% [X,Y] = meshgrid(x,y);
% mult=1;
% X = mult*reshape(X,1,[])';
% Y = mult*reshape(Y,1,[])';
mag = sqrt(U.^2+ V.^2);

quiver(xc, yc, U./mag, V./mag, 'k');
axis equal
hold off

figure(3)

gravedad = (U.^2+V.^2).^0.5;

gravedad = reshape(gravedad,divx,divy);

imagesc(gravedad)

colorbar