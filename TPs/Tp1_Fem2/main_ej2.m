% Main del trabajo practico 1

clc
clear all
close all

%% Datos del problema
E = 6.825e7; %psi
v = 0.3;

R = 10; %in
t = 0.04 ; %in
tita = 90; %grados
phi = 72;

F = 1; %lbf

%% Control
type = 1; %tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

div = 16; %cuantas divisiones en cada  lado

switch type
    case 1 %Mindlin
        T = t;
        dofselem = 6;
        crearK = @crearK_shellMQ4;
        sym_yz = [1 5 6];
        sym_xz = [2 4 6];
        sym_xy = [3 4 5];
    case 2
        T = t*ones(1,4);
        dofselem = 5;
        crearK = @crearK_shell_degeneradoQ4;
        sym_yz = [1 5];
        sym_xz = [2 4];
        sym_xy = [3 4 5];
end
%% Malla

[nodos,elems] = mallador_ej2(phi,R,tita,div,div);

figure(2)
draw_Mesh(elems,nodos, 'NodeLabel',true,'Type','Q4','Color','b')
hold off


nnod = size(nodos,1); %cant nodos
ndof = dofselem*nnod; 
nelem = size(elems,1); %cant elementos

dofs = reshape(1:ndof,dofselem,[])';

%% Calcular direccion z
v3 = direcciones(nodos,elems);
%% Matriz de rigidez

K = zeros(ndof);


for i=1:nelem    
    nodoid = elems(i,:);

    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    v3_el = v3(:,nodoid);
    
    Kel = crearK(nodos(nodoid,:),E,v,T);
    %Kel = crearK(nodos(nodoid,:),E,v,T,v3_el);

    K(dir,dir)=K(dir,dir) + Kel;
end

%% Condiciones de borde
free = true(ndof,1);

borde_BD = abs(nodos(:,1)) <1e-6; %x = 0
borde_AC = abs(nodos(:,2)) <1e-6; %y = 0
borde_AB = 1:(div+1):nnod; %y = L/2
borde_CD = (div+1):(div+1):nnod;

free(dofs(borde_AC,sym_xz)) = false; %sym (mov en x, giro en y z) 6
free(dofs(borde_BD,sym_yz)) = false; %sym (mov y, giro en x y z) 6

%% Cargas

% Definición de las cargas aplicadas
R = zeros(ndof, 1);
A = 1;
B = borde_AB(end);

% eje X
dir = dofs(A,:);
dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas   
dir = dir(1); % Carga solamente en x
R(dir) = R(dir) - F; 

%eje Y
dir = dofs(B,:);
dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas   
dir = dir(2); % Carga solamente en y
R(dir) = R(dir) + F; 



%% Calculos

Kr = K(free,free);

Rr = R(free);

%Calculo movimientos
U = zeros(ndof,1);
U(free) = Kr\Rr;

%% Graficar
x = nodos(:,1);
y = nodos(:,2);
z = nodos(:,3);

escala = 5;
x_deformada = x + escala*U(1:dofselem:ndof);
y_deformada = y + escala*U(2:dofselem:ndof);
z_deformada = z + escala*U(3:dofselem:ndof);
nodos_deformada = [x_deformada y_deformada z_deformada];

figure(3)
draw_Mesh(elems,nodos,'Type','Q4','Color','b')
hold on
draw_Mesh(elems,nodos_deformada,'Type','Q4','Color','k')
hold off

%% Valores del benchmark
Ua = U(dofs(1,1)) %el ultimo nodo en z

%% Tensiones

esfuerzos = zeros(nelem,7); %Nx, Ny, Mx, My, Mxy, Qx, Qy
for i = 1:nelem
    nodoid = elems(i,:);

    dir = dofs(nodoid,:);
    dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

    Coord = nodos(nodoid,:);
    v3_el = v3(:,nodoid);
    
    Uel = U(dir);
    
    esfuerzos(i,:) = stress_shellMQ4(Coord, Uel, E,v,t);
    %esfuerzos(i,:) = stress_shell_degenerado(Coord, Uel, E,v,T,v3_el);
end

% plotear en linea AB
dir = div:div:nelem; %son los ultimos elementos de cada fila
Nx_p = esfuerzos(dir,2); %x' es y en el eje de cordenadas mio
My_p = esfuerzos(dir,3); %y' es x en el eje de cordenadas mio
Qy_p = esfuerzos(dir,6);
% figure(5)
% plot(Nx_p)
% figure(6)
% plot(My_p)
% figure(7)
% plot(Qy_p)
