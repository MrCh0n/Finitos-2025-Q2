clc
clear

auxiliar = [0 0 1];%direccion

De = 40e-3;%m diametro
Di = 30e-3;

%Tubo hueco y acero
E = 210e9;
G = 70e9;
rho = 7850;
Iz = pi*(De^4-Di^4)/64;
Iy = Iz;
A = pi*(De^2-Di^2)/4;
Ip = Iz*2;
y_max = De/2;%m
z_max = De/2;%m
r = De/2;

Nodos = [600 -550 0;
         486 -523 237;
         0 -420 950;
         0 -550 0;
         600 550 0;
         486 523 237;
         50 420 950;
         0 550 0;
         486 -349 237;
         486 349 237;
         288 0 560;
         0 0 560]*1e-3;

elems = [1 2;
         1 9
         2 9
         2 3
         3 7
         3 11
         3 4
         3 12
         4 12
         5 6
         7 8
         5 10
         6 10
         6 7
         7 12
         7 11
         8 12
         9 10
         9 11
         10 11];

nelem = size(elems,1);

div = 4*ones(nelem,1);

[mesh.nodos, mesh.elems.con, mesh.elems.transformada, mesh.elems.inicio] = mesh_1D(elems, Nodos,div);

nnod = size(mesh.nodos,1);
ndof = 6*nnod;%viga 3D
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,6,[])';

mesh.free = true(nnod,6);
empotrados = mesh.elems.transformada([1 4 5 8]);
mesh.free(empotrados,:) = false;

mesh.free = reshape(mesh.free',[],1);

K = zeros(ndof);
M = zeros(ndof);

for i = 1:size(div,1)
    nodoid = mesh.elems.con(mesh.elems.inicio(i), :);

    nodos = mesh.nodos(nodoid,:);

    Kglobal = crearK_viga_3D(nodos,E,G,A,Iy,Iz,Ip, auxiliar);
    Mglobal = masa_viga_3D(nodos,rho,A,Ip,auxiliar);

        %meterlo en la matriz enorme
    cantidad = 0:div(i)-1;
    for j = (mesh.elems.inicio(i) + cantidad)
        dofs1 = mesh.dofs(mesh.elems.con(j,1),:);
        dofs2 = mesh.dofs(mesh.elems.con(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir
        
        K(eleDofs,eleDofs) = K(eleDofs, eleDofs) + Kglobal;
        M(eleDofs,eleDofs) = M(eleDofs,eleDofs) + Mglobal;
        %KelGlobal; 
        %spy(K);
    end 
end


Kr = K(mesh.free,mesh.free);
Mr = M(mesh.free,mesh.free);

[Avec,Aval] = eig(full(Kr),full(Mr));
% Extraer los valores propios y vectores propios
[f_n,idx]= sort(sqrt(diag(Aval)) / (2 * pi)); % Frecuencias naturales en Hz
modos = zeros(sum(mesh.free),ndof); %menos la cantidad de restricciones
modos(:,mesh.free) = Avec(:,idx)'; % Modos de vibración


Frecuencias = f_n(1:6);

%% ej_14
carga = 3;

R = zeros(ndof,1);
switch(carga)%crear R
    case 1%carga destribuida en +Y
        L = norm(Nodos(1,:)-Nodos(2,:)) + norm(Nodos(3,:)-Nodos(2,:)) + norm(Nodos(3,:)-Nodos(4,:));%1-2 2-3 3-4
        q = 40e3/L;%N/m

        pos12 = mesh.elems.inicio(1);%es el elemento 1 originalmente;
        pos23 = mesh.elems.inicio(4);%es el elemento 4 originalmente;
        pos34 = mesh.elems.inicio(7);%es el elemento 7 originalmente;

        pos = pos12:pos12+div(1)-1;
        pos = [pos pos23:pos23+div(4)-1];
        pos = [pos pos34:pos34+div(7)-1];

        pos = mesh.elems.con(pos,:);
        for i = 1:size(pos,1)
            nodoid = pos(i,:);

            nodos = mesh.nodos(nodoid,:);

            L = norm(nodos(2,:)-nodos(1,:));

            dofs1 = mesh.dofs(nodoid(1), :);
            dofs2 = mesh.dofs(nodoid(2), :);

            R(dofs1,:) = R(dofs1,:) + [0 q*L/2 0 -q*L^2/12 0 0]';
            R(dofs2,:) = R(dofs2,:) + [0 q*L/2 0 q*L^2/12 0 0]';
        end 
    case 2%carga en el punto 3;
        pos3 = mesh.elems.transformada(3);
        R(mesh.dofs(pos3,1:3)) = 40e3*[0 1 -1]/sqrt(2);

    case 3%carga destribuida en 3-7
        L =  norm(Nodos(3,:)-Nodos(7,:));
        q = 40e3/0.840;%N/m
        pos = mesh.elems.inicio(5);%es el elemento 5 originalmente;
        pos = pos:pos+div(5)-1;
        pos = mesh.elems.con(pos,:);

        for i = 1:size(pos,1)
            nodoid = pos(i,:);

            nodos = mesh.nodos(nodoid,:);

            L = norm(nodos(2,:)-nodos(1,:));

            dofs1 = mesh.dofs(nodoid(1), :);
            dofs2 = mesh.dofs(nodoid(2), :);

            R(dofs1,:) = R(dofs1,:) + [0 0 -q*L/2 -q*L^2/12 0 0]';
            R(dofs2,:) = R(dofs2,:) + [0 0 -q*L/2 q*L^2/12 0 0]';
        end 
end

Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

Ks = zeros(ndof);

for i = 1:size(div,1)
    nodoid = mesh.elems.con(mesh.elems.inicio(i), :);

    nodos = mesh.nodos(nodoid,:);

    Kglobal = crearK_viga_3D(nodos,E,G,A,Iy,Iz,Ip, auxiliar);

        %meterlo en la matriz enorme
    cantidad = 0:div(i)-1;
    for j = (mesh.elems.inicio(i) + cantidad)
        dofs1 = mesh.dofs(mesh.elems.con(j,1),:);
        dofs2 = mesh.dofs(mesh.elems.con(j,2),:);
        eleDofs = [dofs1 dofs2];%donde tiene que ir

        Uel = U(eleDofs);

        Ks(eleDofs,eleDofs) = Ks(eleDofs, eleDofs) + pandeo_Ks_3D(nodos,Uel, Kglobal,auxiliar);
    end 
end

%Eigenvalue
Ksreducida = Ks(mesh.free,mesh.free);
[Avec_pandeo, Aval_pandeo] = eig(full(Kr), -Ksreducida);

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
modos_final(:,mesh.free) = modos_activos_pandeo'; 

Pandeo = Aval_activos_pandeo(1:6);