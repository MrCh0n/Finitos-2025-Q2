clc 
clear

L = 500;
E = 70000;
A = 100;
alpha = 60; %grados

coord = [0 0;
         L/2 L/2*tand(alpha);
         L 0];

elem = [1 2;
        2 3];

nnod = size(coord,1);
ndof = 2*nnod;
nelem = size(elem,1);

K = zeros(ndof);

dof = [];
for i = 1:nnod
    dof = [dof; i*2-1 i*2];
end

for i = 1:nelem
    V = coord(elem(i,2),:) - coord(elem(i,1),:);
    L = norm(V);
    Kele = A*E/L * [1 -1; -1 1];
    
    ang = V/L;
    T = [ang 0 0;
        0 0 ang];

    Kglobal = T'*Kele*T;

    dir = [dof(elem(i,1),:) dof(elem(i,2),:)];

    K(dir, dir) = K(dir, dir) + Kglobal;
end

Libres = true(nnod,2);
Libres(1,:) = false;
Libres(3,:) = false;

Libres = reshape(Libres',[],1);

Fnodos = zeros(ndof,1);

Fnodos(3) = 1000;

Rreducido = Fnodos(Libres);

Kreducido = K(Libres, Libres);

U = zeros(ndof,1);
U(Libres) = Kreducido\Rreducido;

for i = 1:nelem
    V = coord(elem(i,2),:) - coord(elem(i,1),:);
    L = norm(V);
    B = [-1/L 1/L];
    
    ang = V/L;
    T = [ang 0 0;
        0 0 ang];

    nodo1 = dof(elem(i,1),:);
    nodo2 = dof(elem(i,2),:);

    D = T*[U(nodo1) ; U(nodo2)];

    Stress(i) = E*B*D;
end


a = 50; %amplificar la deformacion

%% Graficar
plot(coord(:,1),coord(:,2),'k*')
hold on
Deformada = coord + a*(reshape(U, 2,nnod))';
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
    Xdef = Deformada(elem(i,:),1);
    Ydef = Deformada(elem(i,:),2);

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