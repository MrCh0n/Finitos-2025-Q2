classdef LST < handle
    properties
        nodes
        elems
        dofs
        free
        bordes
        U
        K
        R
        C
        counts
    end

    methods
        function mesh = LST(bordes,divx,divy,C)
            [mesh.nodes, mesh.elems, mesh.bordes] = mallador_triang_LST(bordes,divx,divy);

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

                Kel = crearK_LST(nodos, mesh.C);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

                mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
            end
        end

        function mesh = armar_R(mesh, i, carga_s, carga_v, type)
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Rel = carga_LST(nodos, carga_s, carga_v, "coordenada", type);

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

        function [escala] = dibujar(mesh,porcien)
            % pasar el porcentaje al que se quiere dibujar la deformada
            % devuelve cuanto se tuvo que escalar
            
            max_U = max(abs(mesh.U));
            porciento = porcien/100;
            escala = mesh.counts.L/max_U*porciento;
            %plot_LST(mesh.nodes,mesh.U,mult);
            
            dofselem = 2;
            ndof = mesh.counts.ndof;
            % sin deformar
            figure(2)
            draw_Mesh(mesh.elems,mesh.nodes, 'NodeLabel',true,'Type','LST','Color','b')
            hold off

            % Deformada
            x = mesh.nodes(:,1);
            y = mesh.nodes(:,2);
 
            x_deformada = x + escala*mesh.U(1:dofselem:ndof);
            y_deformada = y + escala*mesh.U(2:dofselem:ndof);
            nodos_deformada = [x_deformada y_deformada];
            
            figure(3)
            draw_Mesh(mesh.elems,mesh.nodes,'Type','LST','Color','b')
            hold on
            draw_Mesh(mesh.elems,nodos_deformada,'Type','LST','Color','k')
            hold off
        end
    end


end