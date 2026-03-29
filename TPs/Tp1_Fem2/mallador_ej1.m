function [nodos, elems] = mallador_ej1(bordes,divx,divy,R,titamax)
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

function [nodos, elems] = mallador_cuadrado_Q4(bordes,divX,divY)
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
    divX = divX + 1;
    divY = divY + 1;

    disX = linspace(-1,1,divX);
    disY = linspace(-1,1,divY);

    nodos = zeros(divX*divY,2);

    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];

    xi = bordes(:,1);
    yi = bordes(:,2);

    for i = 1:divX
        offset = (i-1)*divY;
        x = disX(i);
        for j = 1:divY
            y = disY(j);
            N = [1 x y x*y]/A;
            nodos(offset+j,:) = [N*xi,N*yi];
        end
    end
    
    nelem = (divY-1)*(divX-1);
    elems = zeros(nelem,4);
    
    padding = divY;%el mismo de la siguiente columna
    for i = 1:divX-1
        offset = (i-1)*divY;
        offset_elems = (i-1)*(divY-1);
        for j = 1:divY-1
            elems(offset_elems+j,:) = [offset+j offset+j+padding offset+j+padding+1 offset+j+1];
        end
    end
end