function [nodos, elems] = mallador_cuadrado_Q8(bordes,divX,divY)
%mallador de elementos Q8 para un cuadrilatero
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

    disX = linspace(-1,1,2*divX-1);%el 2* -1 es para hacer los intermedios
    disY_completa = linspace(-1,1,2*divY-1);
    disY_mitad = linspace(-1,1,divY);

    nnod = divX*(2*divY-1) + (divX-1)*divY;
    nodos = zeros(nnod,2);

    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];

    xi = bordes(:,1);
    yi = bordes(:,2);

    for i = 1:2*divX-1
        x = disX(i);

        if mod(i,2) == 1%los que en y tienen todos los puntos
            offset = (i-1)/2*(3*divY-1);%salta hasta el primero de (x,0)
            for j = 1:2*divY-1
                y = disY_completa(j);
                N = [1 x y x*y]/A;
                nodos(offset+j,:) = [N*xi,N*yi];
            end
        else
            offset = (i-2)/2*(3*divY-1) + (2*divY-1);%salta hasta el primero de (x+dx/2,0) con dx el largo del lado de un cuadrado
            for j = 1:divY
                y = disY_mitad(j);
                N = [1 x y x*y]/A;
                nodos(offset+j,:) = [N*xi,N*yi];
            end
        end
    end
    
    nelem = (divY-1)*(divX-1);
    elems = zeros(nelem,8);
    
    padding_i = 2*divY-1;%los nodos que estan en (0,y) en el Q8
    padding_c = padding_i + divY;%los nodos que estan en (1,y) en el Q8
    for i = 1:divX-1
        offset = (i-1)*(3*divY-1);%salta hasta el primero de (x,0)
        offset_elems = (i-1)*(divY-1);
        for j1 = 1:divY-1
            j = 2*j1-1;%para agarrar solo las esquinas inferiores
            elems(offset_elems+j1,1:4) = [offset+j offset+j+padding_c offset+j+padding_c+2 offset+j+2];%las esquinas
            elems(offset_elems+j1,5:8) = [offset+j1+padding_i offset+j+padding_c+1 offset+j1+padding_i+1 offset+j+1];%los intermedios
        end
    end
end