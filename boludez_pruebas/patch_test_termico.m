clc
clear

k = 1;
t = 1;

L = 2;
W = 2;

n = 30;%cantidad de nodos en x
m = 30;%cantidad de nodos en y

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

% flujo = K*T-Q;
% 
% flujoim = reshape(flujo,m,n);
% 
% imagesc(X,Y,flujoim)
% colorbar

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

x = linspace(0,L,n);
y = linspace(0,W,m);
[X,Y] = meshgrid(x,y);
contour(X,Y,Tim,20)
colorbar

hold on
x = linspace(0,L,divx);
y = linspace(0,W,divy);
[X,Y] = meshgrid(x,y);
X = reshape(X,1,[])';
Y = reshape(Y,1,[])';

quiver(X,Y,U,V);
axis equal
hold off
function [Uel, Vel] = flujo_Q4(coord, Tel, k)
    % isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    
    %iso para extrapolar las tensiones en puntos de Gauss/superonvergentes
    x2 = [-1; 1; 1; -1];
    y2 = [-1; -1; 1; 1];
    A_rs = [ones(4,1) x2 y2 x2.*y2];
    
    r = x1*sqrt(3);%las esquinas de xi eta
    s = y1*sqrt(3);

    %donde quiero la tensiones
    puntos = [-sqrt(1/3) sqrt(1/3)];
     %toda la seccion superior se cambia junta
    
    cant = size(puntos,2);
    flujo_rs = zeros(cant^2,2);% (:,1) xx, 2 yy, 3 xy, 4 zz
    
    %tension en los puntos de Gauss
    for i = 1:cant
        offset = (i-1)*cant;
        for j = 1:cant
            Neta = [0, 1, 0 puntos(j)]/A;
            Nzeta = [0, 0, 1, puntos(i)]/A;

            D = [Neta; Nzeta];

            J = D*coord;

            Bs = J\D;
        
            flujo_rs(offset+j,:) = -k*Bs*Tel; %(qx,qy)
        end
    end
    flujo_rs = flujo_rs([1 2 4 3],:);%queda mal ubicado del for loop puntos(i/j)
   

    Uel = 0;
    Vel = 0;
    
    for j = 1:cant_puntos
        N = [1 r(j) s(j) r(j)*s(j)]/A_rs;
        Uel = Uel + N*flujo_rs(:,1);%sumo todos los flujos en x en los nodos del elemento
        Vel = Vel + N*flujo_rs(:,2);%sumo todos los flujos en y los nodos del elemento
    end

    Uel = Uel/cant_puntos;
    Vel = Vel/cant_puntos;
end



