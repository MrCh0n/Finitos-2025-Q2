clear
clc

Coord = [0 0;
    1 0;
    2 0;
    0 1;
    1 1;
    2 1]; %coordenadas x y de los nodos
%Coord = [0 0; 1 0; 0 1];


Elem = [1 2;
    2 3;
    1 4;
    2 4;
    2 5;
    2 6;
    3 6;
    4 5;
    5 6]; %referir a que nodos estan conectados los elementos
%Elem = [1 2; 2 3; 3 1];


E = 200000; %en MPa

A = 1; % en mm^2 pueden ser matrices (diferente para cada elemento)

%armar la matriz
nnod = size(Coord,1); %cant de nodos
ndof = 2*nnod; %degrees of freedom
nelem = size(Elem,1); %cantidad de elementos
K = zeros(ndof);

for i=1:nelem
    %matriz de rigidez
    V =Coord(Elem(i,2),:)-Coord(Elem(i,1),:); %vector del elemento
    L = norm(V);
    Kel = A*E/L *[1 -1;-1 1];
    
    %giro la matriz
    cos_seno = V/L ; %cosenos del vector
    cos_tita = cos_seno(1);
    sen_tita = cos_seno(2);
    Q = [cos_tita sen_tita 0 0;
        0 0 cos_tita sen_tita];
    KelGlobal = Q'*Kel*Q;
    
    %meterlo en la matriz enorme
    dir = [Elem(i,1)*2-1 Elem(i,1)*2 Elem(i,2)*2-1 Elem(i,2)*2];%donde tiene que ir
    
    K(dir,dir) = K(dir, dir) + KelGlobal;

end


%Condiciones de borde
BC =[1 7 8]; % nodos agarrados (cada nodo tiene 2 elemento)
%BC = [1 5 6];
Libres = 1:ndof;
Libres(BC) = []; %eliminar los que tengan BC

Kreducida = K (Libres, Libres);

%Fuerzas
R =zeros(ndof,1);
R(4) =- 2000;
R(6)=- 2000; %global
% R(4) = -2000;
Rreducida = R(Libres); %reducido

%Calcular movimientos
U = zeros(ndof,1);
U(Libres) = Kreducida\Rreducida; %matriz\vector es inversa de matriz * vector



plot(Coord(:,1),Coord(:,2),'k*')
hold on
Deformada = Coord + (reshape(U, 2,nnod))';
plot(Deformada(:,1),Deformada(:,2),'bo');
hold on
%spy(K)

%graficar elementos 
for i = 1:nelem
%     %reticulado sin deformar
%     Xs = Coord(Elem(i,:),1);
%     Ys = Coord(Elem(i,:),2);
%     plot(Xs,Ys,'k-');
%     hold on
    
    %reticulado deformado
    Xdef = Deformada(Elem(i,:),1);
    Ydef = Deformada(Elem(i,:),2);
    plot(Xdef, Ydef, 'b-');
    hold on
end
hold off

%Calcular fuerzas
F = K*U; 
sum(F); %comprobar que sumatoria de fuerzas y reacciones de 0

%Calcular tensiones (a chequear)
Stress = zeros(ndof,1);
for i=1:nelem
    
    V =Coord(Elem(i,2),:)-Coord(Elem(i,1),:); %vector del elemento
    L = norm(V);
    Bel = [-1/L 1/L];
    
    %giro la matriz
    cos_seno = V/L ; %cosenos del vector
    cos_tita = cos_seno(1);
    sen_tita = cos_seno(2);
    Q = [cos_tita sen_tita 0 0;
        0 0 cos_tita sen_tita];
    dir = [Elem(i,1)*2-1 Elem(i,1)*2 Elem(i,2)*2-1 Elem(i,2)*2];
    DeLocal = Q*U(dir); %elongacion de la barra en coordenadas locales
    
    Stress(i) = E*DeLocal(1)/L; %E*Bel*DeLocal; 
    %Bel*DeLocal es deformacion 
    %multiplicado por E tensiones

end