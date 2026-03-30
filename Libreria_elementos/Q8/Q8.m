classdef Q8 < handle
    properties
        nodes
        elems
        dofs
        free
        U
        K
        R
        C
        counts
    end

    methods
        function mesh = Q8(bordes,divx,divy,C)
            [mesh.nodes, mesh.elems] = mallador_cuadrado_Q8(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodes,1);
            mesh.counts.ndof = mesh.counts.nnod*2;
            mesh.counts.nelem = size(mesh.elems,1);
            mesh.counts.L = norm(bordes(2,:)-bordes(1,:))/divx;

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);

            mesh.K = spalloc(mesh.counts.ndof, mesh.counts.ndof, mesh.counts.ndof*32);
            mesh.C = C;

            mesh.dofs = reshape(1:mesh.counts.ndof,2,[])';
            mesh.free = true(mesh.counts.nnod,2);
        end

        function mesh = armar_K(mesh)
            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Kel = crearK_Q8(nodos, mesh.C);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

                mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
            end
        end

        function mesh = armar_R(mesh, i, carga_s, carga_v, type)
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Rel = carga_Q8(nodos, carga_s, carga_v, "coordenada", type);

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

        function dibujar(mesh)
            max_U = max(mesh.U);
            porciento = 5/100;
            mult = mesh.counts.L/max_U*porciento;
            plot_Q4(mesh.nodes, mesh.U,mult);
        end
    end


end