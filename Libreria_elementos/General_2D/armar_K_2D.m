function [K] = armar_K_2D(nodos, elems, dofs, dofselem, C, crearK)
    nelem = size(elems,1);
    nnz = nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
    
    I = zeros(nnz,1);
    J = zeros(nnz,1);
    V = zeros(nnz,1);
    cont = 1;

    for i = 1:nelem
        nodoid = elems(i,:);

        coord = nodos(nodoid,:);

        Kel = crearK(coord, C);

        dir = reshape(dofs(nodoid,:)',1,[]);

        %mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
        for a = 1:dofselem
            for b = 1:dofselem
                I(cont) = dir(a);
                J(cont) = dir(b);
                V(cont) = Kel(a,b);
                cont = cont+1;
            end%b
        end%a
    end

    K = sparse(I,J,V);
end