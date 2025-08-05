function [] = barras(sistema)
% funcion que crea y grafica reticulado de elementos barra
% sistema es un stuct con 4 campos
%   el primer campo son las coordenadas de los nodos
%       ej: sistema.coord =[0 0; 0 2; 2 0] en [m]
%   el segundo campo son entre que 2 nodos estan los elementos
%       ej: sistema.elem =[1 2; 1 3; 2 3]
%   el tercer campo son las condiciones de borde que hacen U = 0
%       ej: sistema.Bordcond = [1 7 8] esto es Ux = 0 en nodos 1 y 4
%       Uy = 0 en nodo 4
%   el cuarto campo son las fuerzas externas y en que nodo
%       ej: sistema.fuerzas = [4 -2000; 6 -2000] esto es -2000N en Y en
%       los nodos 2 y 3
campos = fieldnames(sistema);

Coord = sistema.(campos{1});
Elem = sistema.(campos{2});
BC = sistema.(campos{3});
Fnodos = sistema.(campos{4});


E = 200e9; %en Pa

A = 50e-6; % en m^2 pueden ser matrices (diferente para cada elemento)

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

Libres = 1:ndof;
Libres(BC) = []; %eliminar los que tengan BC

Kreducida = K (Libres, Libres);

R =zeros(ndof,1);
for i = 1:size(Fnodos,1)
    R(Fnodos(i,1)) =Fnodos(i,2);
end
Rreducida = R(Libres);%reducido

%Calcular movimientos
U = zeros(ndof,1);
U(Libres) = Kreducida\Rreducida; %matriz\vector es inversa de matriz * vector


%% Calcular fuerzas y tensiones
F = K*U; 
sum(F); %comprobar que sumatoria de fuerzas y reacciones de 0

%Calcular tensiones (a chequear)
Stress = zeros(nelem,1);
for i=1:nelem
    
    V =Coord(Elem(i,2),:)-Coord(Elem(i,1),:); %vector del elemento
    L = norm(V);
    Bel = [-1/L 1/L];
    
    %giro la matriz
    cos_seno = V/L ; %cosenos del vector
    cos_tita = cos_seno(1);
    sen_tita = cos_seno(2);
    Rot = [cos_tita sen_tita 0 0; 
            0 0 cos_tita sen_tita];
    nodo1 = [Elem(i,1)*2-1 Elem(i,1)*2];
    nodo2 = [Elem(i,2)*2-1 Elem(i,2)*2];
    
    D_global = [U(nodo1);U(nodo2)];
    
    D_local = Rot*D_global; %elongacion de la barra en coordenadas locales

    Stress(i) = E*Bel*D_local; 

end

a = 50; %amplificar la deformacion

%% Graficar
plot(Coord(:,1),Coord(:,2),'k*')
hold on
Deformada = Coord + a*(reshape(U, 2,nnod))';
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

    if Stress(i)>0.1
        plot(Xdef, Ydef, 'r-');
    elseif Stress(i)<-0.1
        plot(Xdef, Ydef, 'b-');
    else
        plot(Xdef, Ydef, 'g-');
    end
    hold on
end
hold off