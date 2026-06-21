function [nodos, elems] = mallador_ej1_Q8(largo,R,titamax,divY,divTita)
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

    disY_completa = linspace(0,largo,2*divY-1);
    disY_mitad = linspace(0,largo,divY);
    disTita = linspace(0,titamax,2*divTita-1);

    nnod = divTita*(2*divY-1) + (divTita-1)*divY;
    nodos = zeros(nnod,3);

    zmin = R*cosd(titamax);
    for i = 1:2*divTita-1
        tita = disTita(i);
        x = R*sind(tita);
        z = R*cosd(tita) - zmin;
        if mod(i,2) == 1%los que en y tienen todos los puntos
            offset = (i-1)/2*(3*divY-1);%salta hasta el primero de (x,0)
            for j = 1:2*divY-1
                y = disY_completa(j);
                nodos(offset+j,:) = [x,y,z];
            end
        else
            offset = (i-2)/2*(3*divY-1) + (2*divY-1);%salta hasta el primero de (x+dx/2,0) con dx el largo del lado de un cuadrado
            for j = 1:divY
                y = disY_mitad(j);
                nodos(offset+j,:) = [x,y,z];
            end
        end
    end
    
    nelem = (divY-1)*(divTita-1);
    elems = zeros(nelem,8);
    
    padding_i = 2*divY-1;%los nodos que estan en (0,y) en el Q8
    padding_c = padding_i + divY;%los nodos que estan en (1,y) en el Q8
    for i = 1:divTita-1
        offset = (i-1)*(3*divY-1);%salta hasta el primero de (x,0)
        offset_elems = (i-1)*(divY-1);
        for j1 = 1:divY-1
            j = 2*j1-1;%para agarrar solo las esquinas inferiores
            elems(offset_elems+j1,1:4) = [offset+j offset+j+padding_c offset+j+padding_c+2 offset+j+2];%las esquinas
            elems(offset_elems+j1,5:8) = [offset+j1+padding_i offset+j+padding_c+1 offset+j1+padding_i+1 offset+j+1];%los intermedios
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

