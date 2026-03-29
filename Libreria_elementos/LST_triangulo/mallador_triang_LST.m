function [nodos, elems] = mallador_triang_LST(bordes,divX,divY)
    divX = divX + 1;
    divY = divY + 1;

    
    disX = linspace(-1,1,divX);
    disY = linspace(-1,1,divY);

    nodos = zeros(divX*divY,2);

    % isoparametrico Q4
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
    elems = zeros(nelem,6);
    
    padding = divY;%el mismo de la siguiente columna
    for i = 1:divX-1
        offset_elems = (i-1)*(divY-1)*2;
        offset = (i-1)*divY;
        for j = 1:divY-1
            elems(offset_elems+2*j-1,1:3) = [offset+j offset+j+1+padding offset+j+1];%los que se ven para arriba
            elems(offset_elems+2*j,1:3) = [offset+j offset+j+padding offset+j+padding+1];%los que se ven para abajo
        end
    end

    % % Pongo Nodo en borde
    % borde = 28;
    % ids_agua = find(nodos(:,1) == 0);
    % %ordeno para que esten de menor y a mayor y
    % coordenadas_y = nodos(ids_agua, 2);
    % [~, indices_orden] = sort(coordenadas_y);
    % ids_agua = ids_agua(indices_orden);
    % aux = find(nodos(ids_agua,2)>borde);%encuentra los que estan arriba del agua
    % ids_agua(aux) = [];%los saca de los nodos
    % nodos(ids_agua(end),2) = borde;

    %agregar nodos LST
    for i=1:nelem
        nodoid = elems(i,1:3);

        Coord = nodos(nodoid,:); 

        %defino los nodos adicionales
        nodo4 = mean(Coord(1:2,:));
        nodo5 = mean(Coord(2:3,:));
        nodo6 = mean(Coord([1 3],:));
        
        %busco si ya existen

        [ya_existe4, indice4] = ismember(nodo4, nodos, 'rows');
        [ya_existe5, indice5] = ismember(nodo5, nodos, 'rows');
        [ya_existe6, indice6] = ismember(nodo6, nodos, 'rows');

        if ~ya_existe4
            % Si no existe, lo agregamos
            nodos(end+1, :) = nodo4;
            % El nuevo índice es la última fila
            indice4 = size(nodos, 1); 
        end
        if ~ya_existe5
            % Si no existe, lo agregamos
            nodos(end+1, :) = nodo5;
            % El nuevo índice es la última fila
            indice5 = size(nodos, 1); 
        end
        if ~ya_existe6
            % Si no existe, lo agregamos
            nodos(end+1, :) = nodo6;
            % El nuevo índice es la última fila
            indice6 = size(nodos, 1); 
        end
        
        elems(i,4:6) = [indice4 indice5 indice6];


    end

end