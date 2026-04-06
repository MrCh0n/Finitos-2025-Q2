function [nodos, elems, borde] = mallador_triang_CST(bordes,divX,divY)
%mallador de elementos CST para un cuadrilatero
%
%bordes es una matrix de (4,2) de los nodos limitantes al cuadrilatero
%
%divx es la cantidad de elementos en la direccion x
%divy es la cantidad de elementos en la direccion y
%
% Salidas
% nodos es una lista de las coordenadas de los nodos
% elems es una lista de indices a los nodos de los elementos [ (divX)*(divy), 8 ]
% bordes es un struct con los nodos de los 4 lados para hacer las condiciones de borde
%  - borde.lados_12
%  - borde.lados_23
%  - borde.lados_34
%  - borde.lados_41
%
%poner los nodos como:
%    3
%   / \
%  /   \
% /     \
%1 - - - 2 
    divX = divX + 1;
    divY = divY + 1;

    disX = linspace(-1,1,divX);
    disY = linspace(-1,1,divY);
    
    nnod = divX*divY;
    nodos = zeros(nnod,2);

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

    borde.lado_12 = 1:padding:nnod;
    borde.lado_23 = nnod-padding+1:nnod;
    borde.lado_34 = padding:padding:nnod;
    borde.lado_41 = 1:padding;
end