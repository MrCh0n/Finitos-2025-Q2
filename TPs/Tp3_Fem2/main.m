clc
clear
close all

addpath(genpath(pwd+"/Q4"))
addpath(genpath(pwd+"/General_2D"))

%% Datos
%Cant de elementos por eje
divx = 4;
divy = 30;

%Geometria
L = 1;%m
W = 10;%m
t = 1;%m

bordes = [0 0;
          L 0;
          L W;
          0 W];
%Tiempo
%T = 6000;%tiempo de simulacion
T = 6000;
dt = 1;%[s] delta de tiempo
nt = ceil(T/dt);%cuantos pasos da el for loop

%Material
lambda = 40e6;%[Mpa] parametro Lame
mu = 40e6;%[Mpa] parametro Lame
Kdr = lambda + 2/3*mu;%Drain bulk modulus

%Poros
alpha = 0.4;% coeficiente de Biot

M = 250/6*1e6;%[Mpa] modulo de Biot caso 1
%M = 6.06e9;%[Mpa] modulo de Biot caso 2

Pp = 0;%22246;%[Pa] presion uniforme inicial
phi = 0.375;% Porosidad
k = 1.01937e-9;%[m^2] Permeabilidad
mu_f = 1;%[Pa s] Viscosidad del fluido

%Carga
q = 10000;%[N/m] Carga distribuida
%Strain
C = [lambda+2*mu lambda 0;lambda lambda+2*mu 0;0 0 mu];
Czz = lambda*[1,1,0];

%% Mesh
[nodos, elems, bordes] = mallador_cuadrado_Q4(bordes, divx, divy);

nnod = size(nodos,1);
ndof = nnod*2;
nelem = size(elems,1);

dofs = reshape(1:ndof,2,[])';

free = true(nnod, 2);
freeP = true(nnod,1);
%empotrado en U
en_x = [bordes.lado_41 bordes.lado_23];
en_y = [bordes.lado_12]; 

free(en_x,1) = false;
free(en_y,2) = false;

free = reshape(free', 1, []);
%P=0 arriba
arriba = bordes.lado_34;
freeP(arriba) = false;

%% Matrices
[K,Cg,F,Mm,Kp] = crear_Matporos_Q4(nodos,elems,dofs,C,t,alpha,M,k,mu_f,dt);

% M_monio = alfa^2/Kdr * Mm (con M=1)
M_monio = alpha^2/Kdr*Mm*M;

%% R estatico
R = zeros(ndof,1);

Q = q*L/divx;% la fuerza en los nodos
arriba = bordes.lado_34;
arr_y = dofs(arriba,2);

R(arr_y) = R(arr_y) - Q;
R(arr_y([1 end])) = R(arr_y([1 end]))/2; %los bordes tienen la mitad de fuerza
R_est = R;

%% for loop de U-P-U-P
U = zeros(ndof,1);
P = zeros(nnod,1);


Kr = K(free,free);
inv_Kr = inv(Kr);

%presion
matriz_P = Mm/dt + Kp;
inv_P = inv(matriz_P);

%Condicion de no drenado
freeP_nodrenado = true(nnod,1);
[U,P] = iterar(U,P,K,Cg,F,Kp,Mm,M_monio,dt,free,freeP_nodrenado,R_est,4);

Pinicial=reshape(P,divy+1,divx+1);
P0=mean(mean(Pinicial));

%Condicion de drenado
matriz_Pr = matriz_P(freeP,freeP);
inv_Pr = inv(matriz_Pr);

P(~freeP) = 0;
[U,P,tiempo,presion] = iterar(U,P,K,Cg,F,Kp,Mm,M_monio,dt,free,freeP,R_est,nt);
%% Stress
Pel = zeros(nelem,1);
for i = 1:nelem
    nodoid = elems(i,:);
    
    presiones = P(nodoid);
    
    Pel(i) = nodos_a_elem_Q4(presiones);
end
[bruto, n_bruto, n_suave] = stress_2D_poroelasticidad(nodos, elems, dofs, U, Pel, alpha, C, Czz, @stress_Q4_poros, @global_Q4, @elem_a_nodos_Q4);
%% Plot
figure()
plot(tiempo,presion,"Color",'b')
hold on

analitico = zeros(1,nt);

Kd = lambda+2/3*mu;
v = lambda/2/(lambda+mu);
E0 = 2*mu*(1-v)/(1-2*v);
Eedo = lambda+2*mu;
S = 1/M+alpha^2/Eedo;

cv=k/mu_f/S;

z=0;

for i = 1:nt
    tmp = 0;
    t = tiempo(i);
    cte = pi^2*t*cv/W^2/4;
    for k = 1:1000
        tmp = tmp + (-1)^(k-1)/(2*k-1)*cos(pi/2*z/W*(2*k-1))*exp(-(2*k-1)^2*cte);
    end
    analitico(i) = 4/pi*P0*tmp;
end
plot(tiempo,analitico,"Color",'k')

xlabel('Tiempo [s]')
ylabel('Presión de poros [Pa]')
title('Consolidación unidimensional: comparación numérica vs analítica')

legend('Numérico (FEM)', 'Analítico (Terzaghi)', 'Location', 'best')
grid on
hold off
escala = 200;

x = nodos(:,1);
y = nodos(:,2);

x_deformada = x + escala*U(1:2:ndof);
y_deformada = y + escala*U(2:2:ndof);
nodos_deformada = [x_deformada y_deformada];

% figure()
% draw_Mesh(elems, nodos,'Type','Q4','Color','b')
% hold on
% draw_Mesh(elems,nodos_deformada,'Type','Q4','Color','k')
% hold off

%draw_stress(nodos,elems,bruto)

%reshape(P,divy+1,divx+1)
%reshape(U(2:2:end),divy+1,divx+1)
