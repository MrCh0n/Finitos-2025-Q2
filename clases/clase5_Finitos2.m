clc
clear
close all

addpath(genpath(pwd+"/../Libreria_elementos"))
%% Datos problema
E = 1;
t = 1;
v = 0.3;
rho = 1;
%Plain stress
C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

W = 1;
L = 10;

bordes = [0 0;
          L 0;
          L W;
          0 W];

divx = 10;
divy = 1;

%% Q4
cant_puntos = 4;

% creo el mesh y calculo la K
meshQ4 = Q4(bordes,divx,divy,C);

% hacer la carga
carga = zeros(4,4);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];
for i = 1:meshQ4.counts.nelem
   meshQ4.armar_R(i,carga,volumen,true);
end

%rigidez
meshQ4.armar_K;

%masa
meshQ4.armar_M(rho);

% Condiciones de borde
borde_inferior = meshQ4.bordes.lado_12;
meshQ4.cond_borde(borde_inferior(1),3);%apooyo fijo
meshQ4.cond_borde(borde_inferior(end),2);%apoyo simple

% Calculo
meshQ4.calc_U;
meshQ4.calc_frecuencias;

figure()
hold on
color = ['r','g','b','k'];
for i = [1 2 4]
x = meshQ4.vibraciones.modos(i,1:2:end)'+meshQ4.nodos.coordenadas(:,1);
y = meshQ4.vibraciones.modos(i,2:2:end)'+meshQ4.nodos.coordenadas(:,2);
coord = [x y];
draw_Mesh(meshQ4.elems,coord,'Type','Q4','Color',color(i))
end

%% Q8
cant_puntos = 8;

% creo el mesh y calculo la K
meshQ8 = Q8(bordes,divx,divy,C);

% hacer la carga
carga = zeros(4,6);

volumen = [zeros(cant_puntos,1) -ones(cant_puntos,1)];
for i = 1:meshQ8.counts.nelem
   meshQ8.armar_R(i,carga,volumen,true);
end

%rigidez
meshQ8.armar_K;

%masa
meshQ8.armar_M(rho);

% Condiciones de borde
borde_inferior = meshQ8.bordes.lado_12;
meshQ8.cond_borde(borde_inferior(1),3);%apooyo fijo
meshQ8.cond_borde(borde_inferior(end),2);%apoyo simple

% Calculo
meshQ8.calc_U;
meshQ8.calc_frecuencias;

figure()
hold on
color = ['r','g','b','k'];
for i = [1 2 4]
x = meshQ8.vibraciones.modos(i,1:2:end)'+meshQ8.nodos.coordenadas(:,1);
y = meshQ8.vibraciones.modos(i,2:2:end)'+meshQ8.nodos.coordenadas(:,2);
coord = [x y];
draw_Mesh(meshQ8.elems,coord,'Type','Q8','Color',color(i))
end

%% Viga

I = W^3*t/12;% t es el espesor y W la altura
A = W*t;
l = L/divx;%largo de las secciones

nodosx = linspace(0,L,divx+1)';

nnod = size(nodosx,1);
ndof = nnod*2;%y,tita
nelem = divx;
nodosy = zeros(nnod,1);

nodos = [nodosx nodosy];

elems = zeros(nelem,2);
for i = 1:nelem
    elems(i,:) = [i i+1];
end

dofs = reshape(1:ndof,2,[])';
free = true(nnod,2);

free(1,:) = false;
free(end,1) = false;
free = reshape(free', 1, []);

R = zeros(ndof,1);
K = zeros(ndof,ndof);
M = zeros(ndof,ndof);
for i = 1:nelem
    nodosid = elems(i,:);
    
    coord = nodos(nodosid,:);

    Kel = crearK_viga(coord,E,I);
    Mel = masa_viga(coord,rho,A);

    dir = reshape(dofs(nodosid,:)',1,[]);

    K(dir,dir) = K(dir,dir) + Kel;
    M(dir,dir) = M(dir,dir) + Mel;
    R(dir([1 3])) = R(dir([1 3])) - A*l*rho;
end
U = zeros(ndof,1);

Kr = K(free, free);
Rr = R(free);

U(free) = Kr\Rr;

Mr = M(free,free);

[Avec,Aval] = eig(Kr,Mr);
nodos_libres = size(Aval,1);
% Extraer los valores propios y vectores propios
[f_n,idx]= sort(sqrt(diag(Aval)) / (2 * pi)); % Frecuencias naturales en Hz
modos = zeros(nodos_libres,ndof); %menos la cantidad de restricciones
modos(:,free) = Avec(:,idx)'; % Modos de vibración

figure()
hold on
for i = 1:3
x = nodos(:,1);
y = modos(i,1:2:end)'+nodos(:,2);
coord = [x y];
plot(x,y,color(i));
end

%% Resultados
min(U(1:2:end))%viga
min(meshQ4.U(2:2:end))%Q4
min(meshQ8.U(2:2:end))%Q8

fprintf("primeras 3 frecuencias Viga, Q4, Q8")
f_n(1:3)
meshQ4.vibraciones.frecuencias(1:3)
meshQ8.vibraciones.frecuencias(1:3)