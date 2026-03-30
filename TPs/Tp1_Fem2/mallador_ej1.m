function [nodos, elems] = mallador_ej1(largo,R,titamax,divY,divTita)
%mallador de elementos Q4 para un cuadrilatero
%
%bordes es una matrix de (4,2) de los nodos limitantes al cuadrilatero
%
%divx es la cantidad de elementos en la direccion x
%divy es la cantidad de elementos en la direccion y
%
%poner los nodos en bordes como:
% 4 - - - 3
% |       |
% |       |
% |       |
% 1 - - - 2 
    divY = divY + 1;
    divTita = divTita + 1;

    disY = linspace(0,largo,divY);
    disTita = linspace(0,titamax,divTita);

    nodos = zeros(divTita*divY,3);

    zmin = R*cosd(titamax);
    for i = 1:divTita
        tita = disTita(i);
        offset = (i-1)*divY;
        x=R*sind(tita);
        z = R*cosd(tita) - zmin;
        for j = 1:divY
            y = disY(j);
            nodos(offset+j,:) = [x,y,z];
        end
    end
    
    nelem = (divY-1)*(divTita-1);
    elems = zeros(nelem,4);
    
    padding = divY;%el mismo de la siguiente columna
    for i = 1:divTita-1
        offset = (i-1)*divY;
        offset_elems = (i-1)*(divY-1);
        for j = 1:divY-1
            elems(offset_elems+j,:) = [offset+j offset+j+padding offset+j+padding+1 offset+j+1];
        end
    end
end

function [nodos, elems] = mallador_ej1_(largo,divx,divy,R,titamax)
bordes = [0,0;
          R*sind(titamax),0;
          R*sind(titamax),largo;
          0,largo];

[nodos, elems] = mallador_cuadrado_Q4(bordes,divx,divy);
a = max(nodos(:,1));

nnod = size(nodos,1);

%R = 1;
%titamax=40;
zmin = R*cosd(titamax);
for i = 1:nnod
    x = nodos(i,1);
    y = nodos(i,2);

    tita = titamax*(x/a);
    x=R*sind(tita);
    z = R*cosd(tita) - zmin;
    aux_nodos(i,:) = [x y z];
end
nodos = aux_nodos;
end

