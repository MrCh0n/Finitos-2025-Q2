%% Martone Nicolas 19/8/2025
clc
clear

q = 1000;%[N]
M = -1000;%[Nm]

E = 70e9; %[Pa]
b = 0.1; %[m]
I = b^4/12; %seccion cuadeada masisa
A = b^2;

%ubicaciones en el sistema global en [m]
mesh.nodes = [0 0;
              1 0;
              3 0;
              4 1];

%que 2 nodos componen el elemeto
mesh.elem = [1 2;
             2 3;
             3 4];

nnod = size(mesh.nodes,1);
ndof = 3*nnod; %hay 3 dof por elemento
nelem = size(mesh.elem,1);

mesh.dofs = reshape(1:ndof, 3, [])';

%que dofs estan libres
mesh.free = true(nnod,3);
mesh.free(1,:) = false;%empotrado
mesh.free(2,1:2) = false;%fijo
mesh.free(3,1:2) = false;%fijo

mesh.free = reshape(mesh.free', 1, []);%lo paso a vector

%cargas en el sistema
R = zeros(ndof,1);
R(ndof) = M;

for i = 1:2
    nodo1 = mesh.elem(i,1);
    nodo2 = mesh.elem(i,2);

    L = norm(mesh.nodes(nodo2) - mesh.nodes(nodo1));

    dirnod1 = mesh.dofs(nodo1 , :);
    dirnod2 = mesh.dofs(nodo2 , :);

    P = q*L/2;
    M = q*L^2/12;
    
    dir = [dirnod1 dirnod2];

    R(dir) = R(dir) + [0;-P; -M;0; -P; M];
end

%cre la matriz K
K = zeros(ndof);
Kloc = zeros(6);
T = zeros(6);
for i = 1:nelem
    nodo1 = mesh.elem(i,1);
    nodo2 = mesh.elem(i,2);
   
    V = mesh.nodes(nodo2,:) - mesh.nodes(nodo1,:);%direccion del elemento
    L = norm(V);

    Kel_bar = A*E/L *[1 -1;-1 1];%K de barra

    Y1 = 12*E*I/L^3;%K de una viga
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kel_vig = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];

    dir_bar = [1 4];
    dir_vig = [2 3 5 6];

    Kloc(dir_bar, dir_bar) = Kel_bar;
    Kloc(dir_vig, dir_vig) = Kel_vig;

    cs = V/L;%cosenos directores
    c = cs(1);
    s = cs(2);

    Q = [c s 0;%matriz de giro
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    Kglobal = T'*Kloc*T;

    dirnod1 = mesh.dofs(nodo1 , :);
    dirnod2 = mesh.dofs(nodo2 , :);
    
    dir = [dirnod1 dirnod2];
    
    K(dir, dir) = K(dir, dir) + Kglobal;%sumo la K al global
end

Kr = K(mesh.free, mesh.free);
Rr = R(mesh.free);

U = zeros(ndof,1);
U(mesh.free) = Kr\Rr;
%% insiso A
desplasamientos = reshape(U,3,[])'; %desplazamiento por nodos

%% matriz B y tensiones insisos b y d
divisiones = 5;%cuantas divisiones por elemento tengo para tensiones y deformaciones
Bloc = zeros(divisiones,6);
for i = 1:nelem
    nodo1 = mesh.elem(i,1);
    nodo2 = mesh.elem(i,2);

    V = mesh.nodes(nodo2,:) - mesh.nodes(nodo1,:);%vector de direccion
    L = norm(V);

%cosenos directores
    cs = V/L;
    c = cs(1);
    s = cs(2);
  
    %matriz de rotacion
    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    dir = [mesh.dofs(nodo1,:) mesh.dofs(nodo2,:)];
   
    Uloc =  T*U(dir);
    
    %creo la Matriz B
    x = linspace(0,L, divisiones);

    padding = ones(divisiones,1);%para que queden del mismo tamanio

    Baxial = [-padding/L padding/L];%derivada de N
    Bviga = [-(6*L - 12*x)/L^3;
        -(4*L - 6*x)/L^2;
        (6*L - 12*x)/L^3;
        -(2*L - 6*x)/L^2]'; %segunda derivada de N

    
    
    dir_bar = [1 4];
    dir_vig = [2 3 5 6];
    
    y = b/2;
    
    Bloc(1:1:divisiones, dir_bar) = Baxial;
    Bloc(1:divisiones, dir_vig) = -y*Bviga;%tension en la parte superior
    
    Stress = E*Bloc*Uloc;
     maximo = max(abs(Stress));%encuentra el maximo estres
     Xlocal = x(find(abs(abs(Stress)- maximo) < 1));%encuentra en que metro de la viga es ese maximo
end
%% Deformaciones insiso C
Nloc = zeros(2*divisiones,6);


for i = 1:nelem%calcula solo el inclinado, pero puede calculara cualquiera que este en i
    nodo1 = mesh.elem(i,1);
    nodo2 = mesh.elem(i,2);

    V = mesh.nodes(nodo2,:) - mesh.nodes(nodo1,:);%vector de direccion
    L = norm(V);

    %cosenos directores
    cs = V/L;
    c = cs(1);
    s = cs(2);
  
    %matriz de rotacion
    Q = [c s 0;
         -s c 0;
         0 0 1];
    T(1:3, 1:3) = Q;
    T(4:6, 4:6) = Q;

    dir = [mesh.dofs(nodo1,:) mesh.dofs(nodo2,:)];
   
    Uloc =  T*U(dir);
    
        %creo la matrix de N    
    x = linspace(0,L, divisiones);
    
    dir_bar = [1 4];
    dir_vig = [2 3 5 6];
    
    Nel_bar = [L-x; x]'/L;
    Nel_vig = [(L^3 - 3*L*x.^2 + 2*x.^3)/L^3;
               (L^2*x - 2*L*x.^2 + x.^3)/L^2;
               (- 2*x.^3 + 3*L*x.^2)/L^3;
               -(- x.^3 + L*x.^2)/L^2]';

    Nloc(1:2:2*divisiones,dir_bar) = Nel_bar;
    Nloc(2:2:2*divisiones,dir_vig) = Nel_vig;
    
    %crea los desplazamientos u,v locales del elemento

    Udiv = reshape(Nloc*Uloc,2,divisiones)'; %en este punto son las deformaciones locales
    
    %para crear los puntos iniciales de las divisiones
    inicial = mesh.nodes(nodo1,:);
    final = mesh.nodes(nodo2,:);
    deltaL = (final-inicial)/(divisiones-1);
    
    Q = [c s;
         -s c];%para girar al global
    for j = 1:divisiones
        Udiv(j,:) = (Q'*Udiv(j,:)')' + inicial + deltaL*(j-1); %lo giro al global y le sumo su punto inicial
    end
    %ahora esta en el global y en su direccion x,y
end

