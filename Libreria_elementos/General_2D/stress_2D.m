function [bruto, bruto_nodal, suavizado_nodal] = stress_2D(nodos, elems, dofs, U, C, Czz, stress, Global, elem_a_nodo)
    nelem = size(elems,1);
    nnod = size(nodos,1);
    nodelem = size(elems,2);

    bruto = zeros(nelem,7);
    bruto_nodal = zeros(nnod, 7);
    count = zeros(nnod,1);
    
    for i = 1:nelem
        nodoid = elems(i,:);
    
        dir = dofs(nodoid,:);
        dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
    
        Coord = nodos(nodoid,:);
        
        Uel = U(dir);
        
        sigmas = stress(Coord, Uel, C, Czz);
    
        s1 = sigmas(4);
        s2 = sigmas(5);
        s3 = sigmas(6);
    
        bruto(i,1:4) = sigmas([1 2 3 6]);
        bruto(i,5) = max([s1 s2 s3]); %sigma_1
        bruto(i,6) = min([s1 s2 s3]); %sigma_3
        bruto(i,7) = sqrt(((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)/2);%von Mises
    
        stress_el = elem_a_nodo(Coord, bruto(i,:));
    
        bruto_nodal(nodoid,:) = bruto_nodal(nodoid,:) + stress_el;
        count(nodoid,:) = count(nodoid,:) + 1;
    end
    
    Global_matrix = armar_Global(nodelem, nelem, elems, nodos, Global);
    
    suavizado_nodal = Global_matrix\bruto_nodal;
    bruto_nodal = bruto_nodal./count;
end

function M = armar_Global(dofselem, nelem, elems, nodos, Global)
    nnz = nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
    
    I = zeros(nnz,1);
    J = zeros(nnz,1);
    V = zeros(nnz,1);
    cont = 1;

    for i = 1:nelem
        nodoid = elems(i,:);

        coord = nodos(nodoid,:);

        Mel = Global(coord);

        %mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
        for a = 1:dofselem
            for b = 1:dofselem
                I(cont) = nodoid(a);
                J(cont) = nodoid(b);
                V(cont) = Mel(a,b);
                cont = cont+1;
            end%b
        end%a
    end

    M = sparse(I,J,V);
end