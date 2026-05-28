% este codigo ta porque me daba pena borrarlo, pero no sirve de nada

clc
clear
close all

addpath(genpath(pwd+"/Q4"))
addpath(genpath(pwd+"/General_2D"))
%% Datos
%Cant de elementos por eje
divx = 2;
divy = 10;

%Geometria
L = 1;%m
W = 10;%m
t = 1;%m

bordes = [0 0;
          L 0;
          L W;
          0 W];
%Tiempo
tiempo = 1; %s
dt = 1;%[s] delta de tiempo
%Material
lambda = 40e6;%[Mpa] parametro Lame
mu = 40e6;%[Mpa] parametro Lame

%Poros
alpha = 0.4;% coeficiente de Biot
M = 250/6*1e6;%[Mpa] modulo de Biot
Pp = 0;%[Mpa] presion uniforme inicial
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
%empotrado en U
bordes_u = [bordes.lado_41 bordes.lado_23];
free(bordes.lado_12,2) = false;
free(bordes_u,1) = false;

free = reshape(free', 1, []);
%% Matrices
dofselem = 8;

nnz = nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros

I = zeros(nnz,1);
J = zeros(nnz,1);
V = zeros(nnz,1);
cont = 1;

Cg = zeros(ndof, nnod);
F = zeros(nnod, ndof);
Mm = zeros(nnod, nnod);
Kp = zeros(nnod, nnod);
for i = 1:nelem
    nodoid = elems(i,:);

    coord = nodos(nodoid,:);

    Kel = t*crearK_Q4(coord, C);

    Cg_el = crearCg_Q4(coord, t, alpha);

    F_el = crearF_Q4(coord, t, alpha, dt);

    M_el = crearM_Q4(coord, 1/M);

    Kp_el = crearK_Q4_poros(coord, k, mu_f, t);

    dir = reshape(dofs(nodoid,:)',1,[]);

    Cg(dir,nodoid) = Cg(dir,nodoid) + Cg_el;
    F(nodoid,dir) = F(nodoid,dir) + F_el;
    Mm(nodoid, nodoid) = Mm(nodoid, nodoid) + M_el;
    Kp(nodoid, nodoid) = Kp(nodoid, nodoid) + Kp_el;

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
Kr = K(free,free);
inv_Kr = inv(Kr);

%presion
Mm = Mm/dt;
inv_P = inv(Mm+Kp);
%% R estatico
R = zeros(ndof,1);

Q = q*L/divx;% la fuerza en los nodos
arriba = bordes.lado_34;
arr_y = dofs(arriba,2);

R(arr_y) = R(arr_y) - Q;
R(arr_y([1 end])) = R(arr_y([1 end]))/2;%los bordes tienen la mitad de fuerza
R_est = R;

%% Calculo
lista_tiempo = 1:dt:tiempo; 
nt = length(lista_tiempo);

U_tiempo = zeros(ndof,nt);
P_tiempo = zeros(nnod,nt);
ResiduosP = {};

%Condiciones iniciales
U_prev = zeros(ndof,1);
P_prev = ones(nnod,1)*Pp;
for i = lista_tiempo
    % preparar
    U_iter = U_prev;
    P_iter = P_prev;
    Residuos_iter = [];
    for j = 1:30 %cambiar a while dsp        
        %preparar
        U = zeros(ndof,1);
        P = ones(nnod,1)*Pp;

        %% Calculo U
        R = R_est - Cg*P_iter;
        
        Rr = R(free); %esto es F + Q(Pk-1)

        U(free) = inv_Kr*Rr;%si haces un for loop ya esta hecho la inversion de matriz
        
        difU = norm(U-U_iter);
        
        %% Calculo P
        
        P = inv_P*(M*P_prev + F*(U-U_iter));
        
        difP = norm(P-P_iter);
    
        % actualizar
        U_iter = U;
        P_iter = P;
        Residuos_iter(:,end+1) =  [difU,difP];
    end
    %guardar
    U_tiempo(:,i) = U;
    P_tiempo(:,i) = P;
    
    Residuos{i} = Residuos_iter;
end


%% Stress
Pel = zeros(nelem,1);
for i = 1:nelem
    nodoid = elems(i,:);
    
    presiones = P(nodoid);
    
    Pel(i) = nodos_a_elem_Q4(presiones);
end
[bruto, n_bruto, n_suave] = stress_2D_poroelasticidad(nodos, elems, dofs, U, Pel, alpha, C, Czz, @stress_Q4_poros, @global_Q4, @elem_a_nodos_Q4);

%% Plot
escala = 100;

x = nodos(:,1);
y = nodos(:,2);

x_deformada = x + escala*U(1:2:ndof);
y_deformada = y + escala*U(2:2:ndof);
nodos_deformada = [x_deformada y_deformada];

figure()
draw_Mesh(elems, nodos,'Type','Q4','Color','b')
hold on
draw_Mesh(elems,nodos_deformada,'Type','Q4','Color','k')
hold off

%draw_stress(nodos,elems,bruto)


%% funciones
