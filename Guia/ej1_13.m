clc 
clear

E = 200e9;

De = 100e-3;
Di = 0;

I = pi*(De^4-Di^4)/64;
A = pi*(De^2-Di^2)/4;
y_max = De/2;

L = 10;

P = pi^2*E*I/L^2;

nelem = 4;
x = zeros(nelem+1,1);
y = linspace(0,L,nelem+1)';
mesh.nodos = [x y];

primero = 1:nelem;
segundo = 2:nelem+1;
mesh.elems.con = [primero' segundo'];

nnod = size(mesh.nodos,1);
ndof = nnod*3;
nelem = size(mesh.elems.con,1);

mesh.dofs = reshape(1:ndof,3,[])';

mesh.free = true(nnod,3);
empotrados = [1];
mesh.free(empotrados,:) = false;
movil = [nelem+1];
%mesh.free(movil,1) = false;%movil en y

mesh.free = reshape(mesh.free',[],1);


R = zeros(ndof,1);
R(end-1) = -P;

K = zeros(ndof);

for i =1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);
   

    Kel = crearK_viga_barra(nodos,E,A,I);

    dofs1 = mesh.dofs(mesh.elems.con(i,1),:);
    dofs2 = mesh.dofs(mesh.elems.con(i,2),:);
    eleDofs = [dofs1 dofs2];%donde tiene que ir

    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Kel;
end

Kr = K(mesh.free,mesh.free);
Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;

Ks = zeros(ndof);

for i = 1:nelem
    nodoid = mesh.elems.con(i,:);

    nodos = mesh.nodos(nodoid,:);
   

    Kel = crearK_viga_barra(nodos,E,A,I);

    dofs1 = mesh.dofs(mesh.elems.con(i,1),:);
    dofs2 = mesh.dofs(mesh.elems.con(i,2),:);
    eleDofs = [dofs1 dofs2];%donde tiene que ir

    global_loads = Kel*U(eleDofs);

    V = nodos(2,:) - nodos(1,:);
    L = norm(V);

    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0;
         -s c 0;
         0 0 1];

    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    local_loads = T*global_loads;
    P = -local_loads(1);

    Ks(eleDofs,eleDofs) = Ks(eleDofs,eleDofs) + P*pandeo_Ks_viga(nodos);
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

%Pandeo = Aval_activos_pandeo(1:6);

k = Aval_activos_pandeo.^-0.5
