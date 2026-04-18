function [wB,Nx_p,My_p,Qy_p] = ej1(type,div)
%hecho con clases para hacer mas facil el cambio entre Mindlin y Degenerado
%type: tipo de elemento: 1 --> "Mindlin" o 2 --> "Degenerado"

%div: cuantas divisiones en cada  lado
%% Datos del problema
E = 4.32e8; %psi
v = 0;

R = 25; %in
L = 50; %in
t = 0.25 ; %in
tita = 40; %grados

q = -90; %psi

%% Control
switch type
    case 1 %Mindlin
        mesh = Mindlin(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5 6];
        sym_xz = [2 4 6];
    case 2
        mesh = Degenerado(L,R,tita,div,div,E,v,t);
        sym_yz = [1 5];
        sym_xz = [2 4];
end
%% Malla
%ya esta mallado necesito nodos para las BC
nodos = mesh.nodes;
nelem = mesh.counts.nelem;

%% Matriz de rigidez
mesh.armar_K;

%% Condiciones de borde
borde_AC = abs(nodos(:,1)) <1e-6; %x = 0
borde_CD = abs(nodos(:,2)) <1e-6; %y = 0
borde_AB = abs(nodos(:,2) - L/2) <1e-6; %y = L/2

mesh.cond_borde(borde_CD,[1,3,5]);
mesh.cond_borde(borde_AB,sym_xz);
mesh.cond_borde(borde_AC,sym_yz);

%% Cargas

mesh.armar_R(q);

%% Calculos

mesh.calc_U;

esfuerzos = mesh.esfuerzo;

% plotear en linea AB
dir = div:div:nelem; %son los ultimos elementos de cada fila
Nx_p = esfuerzos(dir,2); %x' es y en el eje de cordenadas mio
My_p = esfuerzos(dir,3); %y' es x en el eje de cordenadas mio
Qy_p = esfuerzos(dir,6); 
%% Graficar
%mesh.dibujar;

%% Valores del benchmark
wB = mesh.U(mesh.dofs(mesh.counts.nnod,3));%el ultimo nodo en z
end