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