function [nodos, elems, borde] = mallador_cuadrado_Q4(bordes,divX,divY)
%mallador de elementos Q4 para un cuadrilatero
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
%  - bordes.lados_12
%  - bordes.lados_23
%  - bordes.lados_34
%  - bordes.lados_41
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

    nnod = divX*divY;

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

    borde.lado_12 = 1:padding:nnod;
    borde.lado_23 = nnod-padding+1:nnod;
    borde.lado_34 = padding:padding:nnod;
    borde.lado_41 = 1:padding;
end