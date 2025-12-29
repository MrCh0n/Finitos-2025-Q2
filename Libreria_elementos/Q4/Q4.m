classdef Q4
    properties
        nodes
        elems
        dofs
        free
        U
        K
        R
        counts
    end

    methods
        function mesh = Q4(bordes,divx,divy)
            [mesh.nodes, mesh.elems] = mallador_cuadrado_Q4(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodes,1);
            mesh.counts.ndof = mesh.counts.nnod*2;
            mesh.counts.nelem = size(mesh.elems,1);

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);

            mesh.K = zeros(mesh.counts.ndof, mesh.counts.ndof);

            mesh.dofs = reshape(1:mesh.counts.ndof,2,[])';
            mesh.free = true(mesh.counts.nnod,2);

        end

        function mesh = armar_K(mesh,i,C)
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Kel = crearK_Q4(nodos, C);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

                mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
        end

        function mesh = armar_R(mesh, i, carga_s, carga_v, type)
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Rel = carga_Q4(nodos, carga_s, carga_v, "coordenada", type);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

                mesh.R(dir) = mesh.R(dir) + Rel;
        end

        function mesh = cond_borde(mesh, nodos, restricciion)
            if restricciion == 3
                mesh.free(nodos,:) = false;
            else
                mesh.free(nodos, restricciion) = false;
            end
        end

        function mesh = calc_U(mesh)
            mesh.free = reshape(mesh.free', 1, []);
            Kr = mesh.K(mesh.free, mesh.free);
            Rr = mesh.R(mesh.free);

            mesh.U(mesh.free) = Kr\Rr;
        end
    end


end