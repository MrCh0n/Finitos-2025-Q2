clc
clear

q_max = 300000; %N/m
q_min = 0;

q = (q_max-q_min)/8;

E = 210e9; %Pa
I = 0.5*1e-4; %m^4
A = 0.5*1e-2; %m^2

mesh.nodes = [0 0;
              0 4;
              0 8; 
              10 4;
              10 8];

mesh.elem = [1 2;
             2 3;
             2 4;
             3 5];

nnod = size(mesh.nodes, 1);
ndof = 3*nnod; %3 porque tengo u,v,tita
nelem = size(mesh.elem,1);

mesh.dofs = reshape(1:ndof, 3, [])';

mesh.Libres = true(nnod, 3);
empotrados = [1 4 5];
mesh.Libres(empotrados,:) = false;

mesh.Libres = reshape(mesh.Libres', [], 1);

R = zeros(ndof,1);
for i = 1:2
    nodo1 = mesh.elem(i,1);
    nodo2 = mesh.elem(i,2);
    q_alto = q*(8-mesh.nodes(nodo1, 2));
    q_bajo = q*(8-mesh.nodes(nodo2, 2));
    w = q_alto-q_bajo;
    P1 = (-7*w/20-q_bajo/2)*4;
    P2 = (-3*w/20-q_bajo/2)*4;
    M1 = (-w/20-q_bajo/12)*4^2;
    M2 = (+w/30+q_bajo/12)*4^2;
    
    dir = [mesh.dofs(nodo1,:) mesh.dofs(nodo2,:)];
    R(dir) = R(dir) + [-P1; 0; M1; -P2; 0; M2];
end

T = zeros(6);
K = zeros(ndof);
Kloc = zeros(6);
for i = 1:nelem
    dir_nod = mesh.elem(i,:);

    nodos = mesh.nodes(dir_nod,:);

     V = nodos(2,:) - nodos(1,:);
     L = norm(V);


    Kloc = crearK_viga_barra(nodos, E, A, I);

    cs = V/L;
    c = cs(1);
    s = cs(2);

    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;
    
    Kglob = T'*Kloc*T;

    dir = mesh.dofs(dir_nod,:)';

    K(dir,dir) = K(dir,dir) + Kglob;
end

Rr = R(mesh.Libres);
Kr = K(mesh.Libres, mesh.Libres);

U = zeros(ndof,1);
U(mesh.Libres) = Kr\Rr;


%% Grafico de desplasamientos
divisiones = 20;
mult = 100; %multiplica los desplazamientos

plot_viga_barra(mesh.nodes, mesh.elem,mesh.dofs, U, divisiones,mult)