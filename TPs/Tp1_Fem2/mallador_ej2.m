function [nodos, elems] = mallador_ej2(phimax,R,titamax,divPhi,divTita)
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
    divPhi = divPhi + 1;
    divTita = divTita + 1;

    disPhi = linspace(0,phimax,divPhi);
    disTita = linspace(0,titamax,divTita);

    nodos = zeros(divTita*divPhi,3);

    for i = 1:divTita
        tita = disTita(i);
        offset = (i-1)*divPhi;
        for j = 1:divPhi
            phi = disPhi(j);
            z = R*sind(phi);
            y = R*sind(tita)*cosd(phi);
            x = R*cosd(tita)*cosd(phi);
            nodos(offset+j,:) = [x,y,z];
        end
    end
    
    nelem = (divPhi-1)*(divTita-1);
    elems = zeros(nelem,4);
    
    padding = divPhi;%el mismo de la siguiente columna
    for i = 1:divTita-1
        offset = (i-1)*divPhi;
        offset_elems = (i-1)*(divPhi-1);
        for j = 1:divPhi-1
            elems(offset_elems+j,:) = [offset+j offset+j+padding offset+j+padding+1 offset+j+1];
        end
    end
end

