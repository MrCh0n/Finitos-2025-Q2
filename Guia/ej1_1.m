%1.1 barra distinto area y material
clc 
clear

F = 2000;
E1 = 70e9;
E2 = 110e9;
A1 = 50^2*pi/4/10^6;
A2 = A1*30^2/50^2;

coord = [0 2 4 5 8]/10;

elem = [1 2 E1 A1;
        2 3 E1 A1;
        3 4 E2 A1;
        4 5 E2 A2];

BC = 1;

Fnodos = [2 3*F;
          4 -2*F;
          5 F];

nnod = size(coord , 2);
nelem = size(elem , 1);

K = zeros(nnod);

for i = 1:nelem
    V = coord(elem(i,2)) - coord(elem(i,1));
    E = elem(i,3);
    A = elem(i,4);
    Kel = A*E/V * [1 -1; -1 1];

    dir = [elem(i,1) elem(i,2)];

    K(dir,dir) = K(dir,dir) + Kel;
end

Libres = 1:nnod;
Libres(BC) = [];

R = zeros(nnod,1);
for i = 1:size(Fnodos,1)
    R(Fnodos(i,1)) = Fnodos(i,2);
end

Rreducido = R(Libres);
Kreducido = K(Libres, Libres);

U = zeros(nnod,1);
U(Libres) = Kreducido\Rreducido;

%Calcular tensiones (a chequear)
Stress = zeros(nelem,1);
for i=1:nelem
    
    V =coord(elem(i,2))-coord(elem(i,1)); %vector del elemento
    L = norm(V);
    Bel = [-1/L 1/L];
    
    D_local = [U(elem(i,1));U(elem(i,2))]; %elongacion de la barra en coordenadas locales

    Stress(i) = E*Bel*D_local; 

end

a =1;%amplificar la deformacion

%% Graficar
plot(coord(:),0,'k*')
hold on
Deformada = coord + a*U'
plot(Deformada(:),0,'bo');
%plot([0 0.2], [0 0], 'b-')
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
    Xdef = Deformada([elem(i,1) elem(i,2)]);
    Y = zeros(size(Xdef));

    if Stress(i)>0.1
        plot(Xdef, Y,'r');
    elseif Stress(i)<-0.1
        plot(Xdef, Y,'b');
    else
        plot(Xdef, Y,'g');
    end
    hold on
end
hold off
