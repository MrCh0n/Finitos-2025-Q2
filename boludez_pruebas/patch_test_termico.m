clc
clear

k = 1;
t = 1;

L = 1;
W = 1;

n = 100;
m = 100;

q = 30;%W/m2

divx = n-1;
divy = m-1;
bordes = [0 0; L 0; L W; 0 W];

q_nodo = q*L/divx;

[nodos, elems] = mallador_cuadrado_Q4(bordes, divx, divy);

nnod = size(nodos,1);
ndof = nnod;
nelem = size(elems,1);

izquierda = 1:divy+1;
arriba = divy+1:divy+1:nnod;
abajo = 1:divy+1:nnod;
derecha = nnod-divy:nnod;

K = zeros(ndof,ndof);

dofs = 1:ndof;

Q = zeros(ndof,1);


Q(arriba) = q_nodo;
Q([m nnod]) = Q([m nnod])/2;%esquinas es la mitad

free = true(nnod, 1);

fijos = [izquierda derecha abajo];

free(fijos) = false;

for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = crearK_termico_Q4(coord,k);

    dir = dofs(nodoid);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
        
    K(dir,dir) = K(dir,dir) + Kel;
end


Kr = K(free,free);

Qr = Q(free);

T = zeros(ndof,1);
T(free) = Kr\Qr;


X = nodos(:,1);
Y = nodos(:,2);
Tim = reshape(T,m,n);

imagesc(X,Y,Tim)
colorbar

flujo = K*T-Q;

flujoim = reshape(flujo,m,n);
figure(2)

imagesc(X,Y,flujoim)
colorbar


