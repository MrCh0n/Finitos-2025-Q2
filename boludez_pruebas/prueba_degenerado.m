% Ejercico para cascaras degeneradas
clc
clear all
close all


addpath(genpath(pwd+"/../Libreria_elementos"))

%% Datos del problema
dofselem = 5;
E = 200e9; %Pa
v = 0.3;

a = 1;%2;%1; %m
t = 50e-3; %m
t=[t t t t];

q0 = 1;%10e6; %Pa

%% Ecuacion Analitica
% Cálculo de la deflexión máxima utilizando la ecuación analítica
D = E * t(1)^3 / (12 * (1 - v^2)); % Rigidez de la placa

w0 = -16*q0*a^4/(pi^6*D);

x = linspace(0,a,101);
y = x;

[X,Y ]= meshgrid(x,y);

W = 0;

for m = 1:2:3
    for n = 1:2:3
        const = 1/(m*n*(m^2+n^2)^2);
        W = W + const*w0*sin(m*X*pi/a).*sin(n*Y*pi/a);
    end
end

% Visualización de la deflexión máxima
figure(1)
surf(X, Y, W);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Deflexión (m)');
title('Deflexión Ecuación exacta');
colorbar;
hold off

%% Malla
div_X = 10;
div_Y = 10;

bordes = [0,0;
          0,a;
          a,a;
          a,0];

[nodos, elems] = mallador_cuadrado_Q4(bordes,div_X,div_Y);

% agrego valor en z
z = sqrt(2*a^2-nodos(:,1).^2) -a;
nodos = [nodos z];

aux = elems(:,2);
elems(:,2) = elems(:,4);
elems(:,4) = aux;

% figure(2)
% draw_Mesh(elems,nodos, 'NodeLabel',true,'Type','Q4','Color','b')
% hold off

nnod = size(nodos,1); %cant nodos
ndof = dofselem*nnod;
nelem = size(elems,1); %cant elementos

dofs = reshape(1:ndof,dofselem,[])';

%% Calcular Direccion z
v3_list = direcciones(nodos,elems);

%% Matriz de rigidez
K = spalloc(ndof,ndof, 32*ndof);

for i=1:nelem    
    nodoid = elems(i,:);

    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    

    Kel = crearK_shell_degeneradoQ4(nodos(nodoid,:),E,v,t);
    %Kel = crearK_shellMQ4(nodos(nodoid,:),E,v,t(1));


    K(dir,dir)=K(dir,dir) + Kel;
end

%% Condiciones de borde
free = true(ndof,1);

borde_1 = abs(nodos(:,1)) <1e-6;
borde_2 = abs(nodos(:,2)) <1e-6;
borde_3 = abs(nodos(:,1) - a) <1e-6;
borde_4 = abs(nodos(:,2) - a) <1e-6;

dofs2fix = [dofs(borde_1,:);dofs(borde_2,:);dofs(borde_3,:);dofs(borde_4,:)];

free(dofs2fix) = false;

%% Definición de las cargas aplicadas
R = zeros(ndof, 1);

cargas_v = -q0*ones(4,1);
aux = zeros(dofselem*4,1);
dirs = [3 dofselem-1 dofselem];
dirs = [dirs dirs+dofselem];
dirs = [dirs dirs+dofselem*2];
for i=1:nelem
    nodoid = elems(i,:);
    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    carga = carga_MQ4(nodos(nodoid,[1 2]),cargas_v);
    aux(dirs) = carga;
    R(dir) = R(dir) + aux; % Carga en el borde inferior
end

%% Calculos
Kr = K(free,free);

Rr = R(free);

%Calculo movimientos
U = zeros(ndof,1);
U(free) = Kr\Rr;

%% plotear
X_plot = zeros(nelem,1);
Y_plot = X_plot;
W_plot = X_plot;
for i = 1:nnod
    X_plot(i) = nodos(i,1);
    Y_plot(i) = nodos(i,2);
    W_plot(i) = U(dofs(i,3)); % Asignar la deflexión en Z a W_plot
end

X_plot = reshape(X_plot,div_Y+1,div_X+1);
Y_plot = reshape(Y_plot,div_Y+1,div_X+1);
W_plot = reshape(W_plot,div_Y+1,div_X+1);



figure(4)
surf(X_plot, Y_plot, W_plot);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Deflexión (m)');
title('Deflexión de la Placa de Mindlin');
colorbar;
hold off