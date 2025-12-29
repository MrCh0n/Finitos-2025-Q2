function [nodos, elems] = mallador_triang_CST(bordes,divX,divY)
%mallador de elementos CST para un cuadrilatero
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
    
    nelem = 2*(divY-1)*(divX-1);
    elems = zeros(nelem,3);
    
    padding = divY;%el mismo de la siguiente columna
    for i = 1:divX-1
        offset_elems = (i-1)*(divY-1)*2;
        offset = (i-1)*divY;
        for j = 1:divY-1
            elems(offset_elems+2*j-1,:) = [offset+j offset+j+1+padding offset+j+1];%los que se ven para arriba
            elems(offset_elems+2*j,:) = [offset+j offset+j+padding offset+j+padding+1];%los que se ven para abajo
        end
    end
end