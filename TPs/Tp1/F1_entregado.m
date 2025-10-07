% TP Suspension Formula 1

clear all
clc
close all

flag_no_deformada = false; %por si queres grafico con deformada o no deformada

%Carga = [-6 0 10 0 -1.980 0]*1e3; %Frenado
%Carga = [0 -10 12 -3.3 0 0]*1e3;  %curva 1 I
%Carga = [0 10 12 3.3 0 0]*1e3;  %curva 2 D
Carga = [0 0 15 0 0 0]*1e3;    %Bache

%% Mesh
cantElemxLinea = 5; %subdivisiones de los elementos
divisiones = [cantElemxLinea*ones(1,4) 1   cantElemxLinea   1    cantElemxLinea*ones(1,5)];%primeros 4 son los wishbones, 
             %| wihb      |bar|Roc|Tor|  Direc|Rig    Sin rotulas Tor y Bar son 1
z_hub = 0.33;
z_central = 0.4;

y_par = 0.8; %distancia parche
y_central = 0.17; %distancia al chasis

z_w = 0.1; %distancia wishbones (respecto del centro
z_dir = 0; %distancia direccion
z_sus = 0.05; %distancia a suspension

x_dir = 0.05; %distancia de direccion
x_w = z_w;

L_rocker = 0.03;

Nodos = [0 y_par z_hub+z_w;   %A1
         0 y_par z_hub-z_w;  %B1
         0 y_par z_hub+z_sus; %C
         x_dir y_par z_hub-z_dir; %D
         x_w y_central z_central+z_w; %E
         -x_w y_central z_central+z_w; %F
         x_w y_central z_central-z_w; %G
         -x_w y_central z_central-z_w; %H
         0 y_central z_central-z_w; %I
         0 y_central z_central-z_w+L_rocker; %J1
         x_dir 0.15 z_hub-z_dir; %K
         0 y_par z_hub; %L
         0 y_par z_hub+z_w;   %A2
         0 y_par z_hub-z_w;  %B2
         0 y_central z_central-z_w+L_rocker]; %J2

%coordenadas x y z de los nodos en metros


Elem = [1 5; %A1E
        1 6; %A1F
        2 7; %B1G
        2 8; %B1H
        3 9; %CI
        9 10; %IJ1
        10 15; %J1J2
        4 11; %DK
        13 3; %parche
        12 14; %parche
        3 12; %parche
        4 12]; %parche

% generamos subdivision y transformada es un vector que indica donde fueron
% a parar los nodos originales
[Nodos, Elem, transformada] = crearMesh(Elem, Nodos, divisiones);

nnod = size(Nodos,1); %cant de nodos
ndof = 6*nnod; %degrees of freedom
nelem = size(Elem,1); %cantidad de elementos

Dofs = reshape(1:ndof,6,[])';

% para saber donde comienza cada elemento
inicioElem = zeros(12,1);
inicioElem(1) = 1;
for i = 2:12
    inicioElem(i) = inicioElem(i-1) + divisiones(i-1);
end

%pos son las rotulas: J1,J2 A1,A2 B1,B2

%Rotula en J
pos = [transformada(10);%posicion de los Dofs a copiar
        transformada(15)];%posicion donde se copian los Dofs
Dofs = crearRotula(Dofs, pos, 4);%4 es una rotula en x

%Rotula en A
pos = [transformada(1);
        transformada(13)];

Dofs = crearRotula(Dofs, pos, 4);%4 es una rotula en x

%Rotula en B
pos = [transformada(2);
        transformada(14)]; 

Dofs = crearRotula(Dofs, pos, 4);%4 es una rotula en x

vec = 1:nnod*6; %busca los dofs que sacamos con las rotulas
vec(unique(Dofs)) = [];

for i = flip(vec)%va bajando el valor de los Dofs para que no quede un salto de numeros
    Dofs(Dofs>i) = Dofs(Dofs>i) - 1;
end

ndof = size(unique(Dofs(:)),1);

%% Datos

materiales = struct('nombre', [], 'rho', [], 'E', [],'G', []);
%CFRP
materiales(1).nombre = 'CFRP';
materiales(1).E = 70e9; %en Pa
materiales(1).G = 32e9;
materiales(1).rho = 1600; %kg/m3

%Ti
materiales(2).nombre = 'Ti 6-4';
materiales(2).E = 114e9; %en Pa
materiales(2).G = 42.4e9;
materiales(2).rho = 4500; %kg/m3


%Rigido
materiales(3).nombre = 'Rigido';
materiales(3).E = 1000e15; %en Pa
materiales(3).G = 400e15;
materiales(3).rho = 4500; %kg/m3

%Geometria
%Wishbones
R_max = zeros(8,1);
ri = 15e-3;
re = 30e-3;

R_max(1:4) = re;

sectionData(1).A = pi*(re^2-ri^2);
sectionData(1).Ip = pi/2*(re^4-ri^4);
sectionData(1).Iz = pi/4*(re^4-ri^4);
sectionData(1).Iy = sectionData(1).Iz;
sectionData(1).Rmax =re;

%Barra Suspension
re = 15e-3;
R_max(5) = re;
sectionData(2).A = pi*re^2;
sectionData(2).Ip = 0;
sectionData(2).Iz = 0;
sectionData(2).Iy = 0;
sectionData(2).Rmax =re;

% Rocker
re = 15e-3;
R_max(6) = re;
sectionData(3).A = pi*re^2;
sectionData(3).Ip  = pi/2*re^4;
sectionData(3).Iz = pi/4*re^4;
sectionData(3).Iy = sectionData(3).Iz;
sectionData(3).Rmax =re;

%Barra Direccion
re = 15e-3;
R_max(8) = re;
sectionData(4).A = pi*re^2;
sectionData(4).Ip  = pi/2*re^4;
sectionData(4).Iz = pi/4*re^4;
sectionData(4).Iy = sectionData(4).Iz;
sectionData(4).Rmax =re;

% Resorte Torsion
re = 15e-3;
R_max(7) = re;
Ipt = pi/2*re^4;
Gt= materiales(2).G;
rhot = materiales(2).rho;
Lt = 0.2;
Kt = Gt*Ipt/Lt;
Mt = pi*re^2*Lt*rhot/2;


%Fuerzas
R = zeros(ndof,1);
posL = transformada(12);
R(Dofs(posL,:)) = Carga;


%BC      Solo giro en X para E F G H                    empotramientos para K y J2
BC = [reshape(Dofs(transformada(5:8),[1:3 5:6])', 1, []) Dofs(transformada(11), :) Dofs(transformada(15), :)];

%% Armar la matriz

K = spalloc(ndof,ndof, 64*ndof);
M = spalloc(ndof,ndof, 64*ndof);


%      Wishbones---Rocker----Rigidos----
material = [1 1 1 1 2 2 0 2 3 3 3 3];
geometria = [1 1 1 1 2 3 0 4 1 1 1 1];
for i=[1:4 6 8:12]
    nodo1ID = Elem(inicioElem(i),1);%agarra el primer nodo y su siguiente en los subelementos
    nodo2ID = Elem(inicioElem(i),2);

    nodo1 = Nodos(nodo1ID,:);
    nodo2 = Nodos(nodo2ID,:);
    L = norm(nodo2-nodo1);

    Kel = zeros(12);
    Mel = zeros(12);
    
    
    E=materiales(material(i)).E;
    G=materiales(material(i)).G;
    rho=materiales(material(i)).rho;

    A = sectionData(geometria(i)).A;
    Iz= sectionData(geometria(i)).Iz;
    Iy= sectionData(geometria(i)).Iy;
    Ip= sectionData(geometria(i)).Ip;


    Kel([1 7],[1 7]) = A*E/L*[1 -1;-1 1];
    Kel([2 6 8 12],[2 6 8 12]) = rigidez_viga(L,E,Iz);
    Kel([3 5 9 11],[3 5 9 11]) = rigidez_viga(L,E,Iy);
    Kel([4 10],[4 10]) = G*Ip/L*[1 -1;-1 1];

    Mel([1 7],[1 7]) = rho*A*L/2*[1 0; 0 1];
    Mel([2 6 8 12],[2 6 8 12]) = rho*A*L/24 * diag([12 L^2 12 L^2]);
    Mel([3 5 9 11],[3 5 9 11]) = rho*A*L/24 * diag([12 L^2 12 L^2]);
    Mel([4 10],[4 10]) = rho*Ip*L/2*[1 0; 0 1];

    %giro la matriz
    dir1 = (nodo2-nodo1)/L;
    auxiliarVector = [0 1 0];
    if i == 8 %esta alineado con Y
        auxiliarVector = [1 0 0];
    end

    Q = crearQ(dir1, auxiliarVector);
    
    KelGlobal = Q'*Kel*Q;
    MelGlobal = Q'*Mel*Q;

    %meterlo en la matriz enorme
    cantidad = 0:divisiones(i)-1;
    for j = (inicioElem(i) + cantidad)
        dofs1 = Dofs(Elem(j,1),:);
        dofs2 = Dofs(Elem(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        
        K(eleDofs,eleDofs) = K(eleDofs, eleDofs) + KelGlobal;
        M(eleDofs,eleDofs) = M(eleDofs,eleDofs) + MelGlobal;
        %KelGlobal; 
        %spy(K);
    end

end


% Barra de Suspension------------------------------------------
i = 5;
geo = 2;

nodo1ID = Elem(inicioElem(i),1);
nodo2ID = Elem(inicioElem(i),2);

nodo1 = Nodos(nodo1ID,:);
nodo2 = Nodos(nodo2ID,:);

L = norm(nodo2-nodo1);

Kel = zeros(12);
Mel = zeros(12);

E = materiales(material(i)).E;
rho = materiales(material(i)).rho;

A = sectionData(geometria(i)).A;

Kel([1 7],[1 7]) = A*E/L*[1 -1;-1 1];
Mel([1 7],[1 7]) = rho*A*L/2*[1 0; 0 1];

%giro la matriz
dir1 = (nodo2-nodo1)/L;
auxiliarVector = [0 1 0];

Q = crearQ(dir1, auxiliarVector);

KelGlobal = Q'*Kel*Q;
MelGlobal = Q'*Mel*Q;

j = inicioElem(i);
dofs1 = Dofs(Elem(j,1),:);
dofs2 = Dofs(Elem(j,2),:);
eleDofs = [dofs1 dofs2];

K(eleDofs,eleDofs) = K(eleDofs, eleDofs) + KelGlobal;
M(eleDofs,eleDofs) = M(eleDofs,eleDofs) + MelGlobal;

% Resorte de torsion------------------------------------------
i=7;

Kel = zeros(12);
Mel = zeros(12);

Kel([4 10],[4 10]) = Kt*[1 -1;-1 1];
Mel([4 10],[4 10]) = Mt*[1 0; 0 1];


%no hace falta girarlo
KelGlobal = Kel;
MelGlobal = Mel;

j = inicioElem(i);
dofs1 = Dofs(Elem(j,1),:);
dofs2 = Dofs(Elem(j,2),:);
eleDofs = [dofs1 dofs2];

K(eleDofs,eleDofs) = K(eleDofs, eleDofs) + KelGlobal;
M(eleDofs,eleDofs) = M(eleDofs,eleDofs) + MelGlobal;

%% Calcular Dezplacamientos
%Condiciones de borde
Libres = 1:ndof;
Libres(BC) = []; %eliminar los que tengan BC


Kreducida = K (Libres, Libres);
Rreducida = R(Libres);

U = zeros(ndof,1);

U(Libres) = Kreducida\Rreducida;

R = K*U;

%% Tensiones
plot3(Nodos(:,1),Nodos(:,2),Nodos(:,3),'ko')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on

U_unicos = U;

U_unicos(Dofs(transformada(13:15),4)) = [];%pasamos a los nodos iniciales, las rotaciones quedan al final, 
% el 4 es porque todas estan en TitaX

Nodos_unicos = Nodos;
Nodos_unicos(transformada(13:15),:) = [];

mov = reshape(U_unicos(1:ndof-3),6,[])';%le sacamos las rotaciones A2 B2 J2

Deformada = Nodos_unicos + 40*mov(:,1:3);

plot3(Deformada(:,1),Deformada(:,2),Deformada(:,3),'b*')

div = 5;%divisiones por subelemento
stress = zeros(nelem,4);
%vigas 3D
for i=[1:4 6 8]
    nodo1ID = Elem(inicioElem(i),1);%agarra el primer nodo y su siguiente en los subelementos
    nodo2ID = Elem(inicioElem(i),2);

    nodo1 = Nodos(nodo1ID,:);
    nodo2 = Nodos(nodo2ID,:);
    L = norm(nodo2-nodo1);
    
    %giro la matriz
    dir1 = (nodo2-nodo1)/L;
    auxiliarVector = [0 1 0];
    if i == 8 %esta alineado con Y
        auxiliarVector = [1 0 0];
    end

    Q = crearQ(dir1, auxiliarVector);
    
    %flexion
    x = linspace(0,L,div);
    Bflexion = [-6/L^2+12*x/L^3
                -4/L+6*x/L^2
                 6/L^2-12*x/L^3
                -2/L+6*x/L^2];

    E_elem = materiales(material(i)).E;
    G_elem = materiales(material(i)).G;
    r = sectionData(geometria(i)).Rmax;%radio de los cilindros
    y_max = r;

    cantidad = 0:divisiones(i)-1;
    for j = (inicioElem(i) + cantidad)
        dofs1 = Dofs(Elem(j,1),:);
        dofs2 = Dofs(Elem(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir

        U_local = Q*U(eleDofs);

        stressAx = E_elem*(U_local(7)-U_local(1))/L;

        strain1 = Bflexion'*U_local([2 6 8 12]);   %1 es con v1, titaz1, v2, titaz2
        stress1 = E_elem*y_max*strain1;
    
        strain2 = Bflexion'*U_local([3 5 9 11]);   %2 es con w1, titay1, w2, titay2
        stress2 = E_elem*y_max*strain2;

        stressFlex = max(sqrt(stress1.^2 + stress2.^2));

        stressxx = max(abs([stressFlex+stressAx stressFlex-stressAx]));

        stressxy = G_elem*r*(U_local(10)-U_local(4))/L;
        
        vonMinses = sqrt(stressxx^2 + 3*stressxy^2);

        stress(j, 1:5) = [vonMinses, stressAx, stressFlex,stressxx, stressxy]/1e6;%los paso a Mpa
    end
end

%barra
i=5;
cantidad = 0:divisiones(i)-1;

nodo1ID = Elem(inicioElem(i),1);
nodo2ID = Elem(inicioElem(i),2);

nodo1 = Nodos(nodo1ID,:);
nodo2 = Nodos(nodo2ID,:);

L = norm(nodo2-nodo1);

%giro la matriz
dir1 = (nodo2-nodo1)/L;
auxiliarVector = [0 1 0];

Q = crearQ(dir1, auxiliarVector);

for j = (inicioElem(i) + cantidad)
    dofs1 = Dofs(Elem(j,1),:);
    dofs2 = Dofs(Elem(j,2),:);
    eleDofs = [dofs1 dofs2];%donde tiene que ir

    U_local = Q*U(eleDofs);

    E = materiales(material(i)).E;

    stressAx = E*(U_local(7)-U_local(1))/L;

    stress(j,1:5) = [abs(stressAx), stressAx, 0, 0, 0]/1e6;
end

%ploteo de colores

for i = 1:12
    cantidad = 0:divisiones(i)-1;
    for j = (inicioElem(i) + cantidad)
        nodo1ID = Elem(j,1);
        nodo2ID = Elem(j,2);
       
        if flag_no_deformada
            Xplt = Nodos([nodo1ID, nodo2ID],1);
            Yplt = Nodos([nodo1ID, nodo2ID],2);
            Zplt = Nodos([nodo1ID, nodo2ID],3);
        else
            % bajarle el indice a todos los indices que esten por arriba de
            % las 3 rotulas que sacamos
            
            if nodo1ID==transformada(13)
                nodo1IDcorregido = transformada(1);
            elseif nodo1ID==transformada(14)
                nodo1IDcorregido = transformada(2);
            elseif nodo1ID==transformada(15)
                nodo1IDcorregido = transformada(10);
            else
                nodo1IDcorregido = nodo1ID - sum(nodo1ID>transformada(13:15));
            end
            if nodo2ID==transformada(13)
                nodo2IDcorregido = transformada(1);
            elseif nodo2ID==transformada(14)
                nodo2IDcorregido = transformada(2);
            elseif nodo2ID==transformada(15)
                nodo2IDcorregido = transformada(10);
            else
                nodo2IDcorregido = nodo2ID - sum(nodo2ID>transformada(13:15));
            end
           
            
            Xplt = Deformada([nodo1IDcorregido, nodo2IDcorregido],1);
            Yplt = Deformada([nodo1IDcorregido, nodo2IDcorregido],2);
            Zplt = Deformada([nodo1IDcorregido, nodo2IDcorregido],3);
        end

        if stress(j,2)>0.1%plotea el axial
            plot3(Xplt, Yplt, Zplt, 'r-');
        elseif stress(j,2)<-0.1
            plot3(Xplt, Yplt, Zplt, 'b-');
        else
            plot3(Xplt, Yplt, Zplt, 'g-');
        end
    end
end

title('Deformación de la Estructura (x40)')
h1 = plot3(nan, nan, nan, 'r', 'LineWidth',1.5);
h2 = plot3(nan, nan, nan, 'g', 'LineWidth',1.5);
h3 = plot3(nan, nan, nan, 'b', 'LineWidth',1.5);

legend([h1 h2 h3], {'Tracción','Sin carga axial','Compresión'},'Location','northeast')

stress = stress;%el primer 0 es el resorte los ultimos 4 son los rigidos


%% Frecuencias Naturales

Mreducida = M (Libres, Libres);

[Avec,Aval] = eig(full(Kreducida),full(Mreducida));
% Extraer los valores propios y vectores propios
[f_n,idx]= sort(sqrt(diag(Aval)) / (2 * pi)); % Frecuencias naturales en Hz
modos = zeros(ndof-size(BC,2),ndof); %menos la cantidad de restricciones
modos(:,Libres) = Avec(:,idx)'; % Modos de vibración


Frecuencias = f_n(1:6);

% Graficar frecuencias
figure
plot3(Nodos(:,1),Nodos(:,2),Nodos(:,3),'ko')
title('Modos de Vibraciones')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on
c = ['b*';'r*';'g*';'y*';'c*'];
for ifrecuencias = 1:3

    Ufrecuencias = modos(ifrecuencias,:);
    Ufrecuencias(Dofs(transformada(13:15),4)) = []; %eliminar doble giros
    UfrecuenciasMatriz = reshape(Ufrecuencias(1:ndof-3), 6, [])';
    escala = 0.05./max(abs(UfrecuenciasMatriz(:, 1:3)), [], 'all');
    Deformadafrecuencias = Nodos_unicos + escala*UfrecuenciasMatriz(:,1:3);
    Xfrecuencias = Deformadafrecuencias(:,1);
    Yfrecuencias = Deformadafrecuencias(:,2);
    Zfrecuencias = Deformadafrecuencias(:,3);
    
    plot3(Xfrecuencias,Yfrecuencias,Zfrecuencias,c(ifrecuencias,:))
    
end


for i=1:3
    h(i) = plot3(nan, nan, nan, c(i,:), 'LineWidth',1.5);
    str{i} = sprintf('Modo %d Vibraciones: %.3f Hz',i,Frecuencias(i));
end

legend(h, str,'Location','northeast')
hold off


%% Pandeo

%Conseguir Cargas en los elementos
elementLoads = zeros(nelem,12);
% Solo necesitamos para Wishbones, Rocker y Direccion
for i=[1:4 6 8]
    nodo1ID = Elem(inicioElem(i),1);%agarra el primer nodo y su siguiente en los subelementos
    nodo2ID = Elem(inicioElem(i),2);

    nodo1 = Nodos(nodo1ID,:);
    nodo2 = Nodos(nodo2ID,:);
    L = norm(nodo2-nodo1);

    E = materiales(material(i)).E;
    G = materiales(material(i)).G;

    A = sectionData(geometria(i)).A;
    Iz= sectionData(geometria(i)).Iz;
    Iy= sectionData(geometria(i)).Iy;
    Ip= sectionData(geometria(i)).Ip;

    Kel = zeros(12);
    Kel([1 7],[1 7]) = A*E/L*[1 -1;-1 1];
    Kel([2 6 8 12],[2 6 8 12]) = rigidez_viga(L,E,Iz);
    Kel([3 5 9 11],[3 5 9 11]) = rigidez_viga(L,E,Iy);
    Kel([4 10],[4 10]) = G*Ip/L*[1 -1;-1 1];
    
    %giro la matriz
    dir1 = (nodo2-nodo1)/L;
    auxiliarVector = [0 1 0];
    if i == 8 %esta alineado con Y
        auxiliarVector = [1 0 0];
    end

    Q = crearQ(dir1, auxiliarVector);

    KelGlobal = Q'*Kel*Q;

    cantidad = 0:divisiones(i)-1;
    for j = (inicioElem(i) + cantidad)
        dofs1 = Dofs(Elem(j,1),:);
        dofs2 = Dofs(Elem(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        elementLoads(j, :) = (KelGlobal*U(eleDofs))';
    end
end

%Armar Ks
Ks = zeros(ndof);

for i=[1:4 6 8]
    nodo1ID = Elem(inicioElem(i),1);%agarra el primer nodo y su siguiente en los subelementos
    nodo2ID = Elem(inicioElem(i),2);

    nodo1 = Nodos(nodo1ID,:);
    nodo2 = Nodos(nodo2ID,:);
    L = norm(nodo2-nodo1);

    Ksel = zeros(12);
    Ksel([2 6 8 12],[2 6 8 12]) = 1/(30*L)*[36 3*L -36 3*L;
                                          3*L 4*L^2 -3*L -L^2;
                                          -36 -3*L 36 -3*L;
                                          3*L -L^2 -3*L 4*L^2];
    Ksel([3 5 9 11],[3 5 9 11]) = 1/(30*L)*[36 3*L -36 3*L;
                                          3*L 4*L^2 -3*L -L^2;
                                          -36 -3*L 36 -3*L;
                                          3*L -L^2 -3*L 4*L^2];
    
    %giro la matriz
    dir1 = (nodo2-nodo1)/L;
    auxiliarVector = [0 1 0];
    if i == 8 %esta alineado con Y
        auxiliarVector = [1 0 0];
    end

    Q = crearQ(dir1, auxiliarVector);

    KselGlobal = Q'*Ksel*Q; %esto es sin P

    cantidad = 0:divisiones(i)-1;
    for j = (inicioElem(i) + cantidad)
        dofs1 = Dofs(Elem(j,1),:);
        dofs2 = Dofs(Elem(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        elementalLocalLoad = (Q*elementLoads(j, :)')'; % paso a locales para ver la carga axial
        P = -elementalLocalLoad(1);
        Ks(eleDofs, eleDofs) = Ks(eleDofs, eleDofs) + KselGlobal*P;
    end
end

%Eigenvalue
Ksreducida = Ks(Libres,Libres);
[Avec_pandeo, Aval_pandeo] = eig(full(Kreducida), -Ksreducida);

% Eliminamos infinitos
id_activos_pandeo=find(~isinf(diag(Aval_pandeo)));
Avec_pandeo_activo = Avec_pandeo(:,id_activos_pandeo);

% Ordenar y filtrar
[Aval_activos_pandeo, Index] = sort(diag(Aval_pandeo(id_activos_pandeo,id_activos_pandeo)));
modos_activos_pandeo = Avec_pandeo_activo(:, Index);



%saco negativos
I_positivos = Aval_activos_pandeo>0;
Aval_activos_pandeo = Aval_activos_pandeo(I_positivos);
modos_activos_pandeo = modos_activos_pandeo(:, I_positivos);


modos_final = zeros(size(Aval_activos_pandeo,1),ndof); 
modos_final(:,Libres) = modos_activos_pandeo'; 

Pandeo = Aval_activos_pandeo(1:6);

% Graficar pandeo
figure
plot3(Nodos(:,1),Nodos(:,2),Nodos(:,3),'ko')
title('Modos de Pandeo')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on
c = ['b*';'r*';'g*';'y*';'c*'];

for ipandeo = 1:3
    
    Upandeo = modos_final(ipandeo,:);
    Upandeo(Dofs(transformada(13:15),4)) = []; %eliminar doble giros
    UpandeoMatriz = reshape(Upandeo(1:ndof-3), 6, [])';
    escala = 0.05./max(abs(UpandeoMatriz(:, 1:3)), [], 'all');
    Deformadapandeo = Nodos_unicos + escala*UpandeoMatriz(:,1:3);
    Xpandeo = Deformadapandeo(:,1);
    Ypandeo = Deformadapandeo(:,2);
    Zpandeo = Deformadapandeo(:,3);
    
    plot3(Xpandeo,Ypandeo,Zpandeo,c(ipandeo,:))
    
end

for i=1:3
    h(i) = plot3(nan, nan, nan, c(i,:), 'LineWidth',1.5);
    str{i} = sprintf('Factor de Pandeo %d: %.3f',i,Pandeo(i));
end

legend(h, str,'Location','northeast')
hold off

 
%% Imprimir en pantalla
dofs = Dofs(Elem(inicioElem(10),1),:);
movL = U(dofs);
fprintf('movimiento del chasis en Z [mm]\n%f\n', movL(3)*1e3)
fprintf('\n     Carga (kN)      \n')
fprintf('Fx\t Fy\t Fz\t Mx\t My\t Mz\n')
Carga = Carga/1000;
fprintf('%.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n',Carga(1),Carga(2),Carga(3),Carga(4),Carga(5),Carga(6))

fprintf('\n          \t\tTensiones (MPa)      \n')
fprintf('Elemento\t Tipo\t Ubicacion\t SigmaAx\t SigmaFl\t Sigmaxx\t Sigmaxy\t Von Mises\n')
Elementos = ['WB AE ';'WB AF ';'WB BG ';'WB BH '; 'Suspen'; 'Rocker';'xxxxxx'; 'Direcc'];
Tipo = ['Viga ';'Viga ';'Viga ';'Viga ';'Barra';'Viga ';'xxxxx';'Viga '];
for iImprimir=[1:6 8]
    
    id_stress = sum(divisiones(1:iImprimir-1))+1:sum(divisiones(1:iImprimir)); %conseguir indices que estan las tensiones del elementos
    Tensiones = stress(id_stress,:);
    [~,I] = max(stress(id_stress,1)); %encontrar ubicacion dentro del elemento con maximo
    TensionMax = Tensiones(I,:);

    fprintf('%s \t\t %s \t %d/%d\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\n',Elementos(iImprimir,:),Tipo(iImprimir,:),I,divisiones(iImprimir),TensionMax(2),TensionMax(3),TensionMax(4),TensionMax(5),TensionMax(1))
end

fprintf('\n  Frecuencias Naturales (Hz)      \n')
for ifrec=1:size(Frecuencias,1)
    fprintf('%d\t %4.4f\n',ifrec,Frecuencias(ifrec))
end


fprintf('\n  Pandeo (lambda_cr)      \n')
for ipandeo=1:size(Pandeo,1)
    fprintf('%d\t %4.4f\n',ipandeo,Pandeo(ipandeo))
end


function [Q] = crearQ(dir1, auxiliar)
    dir3 = cross(auxiliar, dir1)/norm(cross(auxiliar, dir1));
    dir2 = cross(dir3, dir1);
    lambda = [dir1; dir2; dir3];
    assert(abs(det(lambda))-1 < 1e-6);
    
    Q = blkdiag(lambda, lambda, lambda, lambda);
end

function [K] = rigidez_viga(L, E, I)
A = E*I/L;
Y1 = A*12/(L^2);
Y2 = A*6/L;
Y3 = A*4;
Y4 = A * 2; 

K =[Y1 Y2 -Y1 Y2;
    Y2 Y3 -Y2 Y4;
    -Y1 -Y2 Y1 -Y2;
    Y2 Y4 -Y2 Y3];
end

function [nodos, elems, transformada] = crearMesh(Elems, Nodes, Divisions)
   
nelem = size(Elems,1);
nnode = size(Nodes,1);

pointsNodeID = zeros(nnode, 1);
Puntocreado = false(nnode,1);

elems = [];
nodos = [];

contador = 1;

for i = 1:nelem
    inicio = Nodes(Elems(i,1), :);
    final = Nodes(Elems(i,2), :);
    
    subnodos = zeros(Divisions(i)+1,3);

    subnodos(:,1) = linspace(inicio(:,1), final(:,1), Divisions(i) +1);
    subnodos(:,2) = linspace(inicio(:,2), final(:,2), Divisions(i) +1);
    subnodos(:,3) = linspace(inicio(:,3), final(:,3), Divisions(i) +1);

    ultimo = contador + Divisions(i);

    subelems = [(contador:ultimo-1)' (contador+1:ultimo)'];

    if Puntocreado(Elems(i,1))
        subnodos(1,:) = [];
        subelems = subelems-1;
        subelems(1,1) = pointsNodeID(Elems(i,1));
    else
        pointsNodeID(Elems(i,1)) = subelems(1,1);
    end

     if Puntocreado(Elems(i,2))
        subnodos(end,:) = [];
        subelems(end,2) = pointsNodeID(Elems(i,2));
    else
        pointsNodeID(Elems(i,2)) = subelems(end,2);
     end
    
    contador = max(max(subelems)) + 1;
    Puntocreado([Elems(i,1) Elems(i,2)]) = true;

    nodos = [nodos; subnodos];
    elems = [elems; subelems];

    transformada(Elems(i,1)) = [subelems(1,1)];
    transformada(Elems(i,2)) = [subelems(end,2)];
end

end

function [Dofs] = crearRotula(Dofs, pos, indice) 
    vec = 1:6;
    vec(indice) = [];
    Dofs(pos(2),vec) = Dofs(pos(1),vec); %solo tita del indice
end
