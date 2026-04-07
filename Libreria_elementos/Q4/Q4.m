classdef Q4 < handle
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
        function mesh = Q4(bordes,divx,divy,C)
            [mesh.nodes, mesh.elems, mesh.bordes] = mallador_cuadrado_Q4(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodes,1);
            mesh.counts.ndof = mesh.counts.nnod*2;
            mesh.counts.nelem = size(mesh.elems,1);
            mesh.counts.L = norm(bordes(2,:)-bordes(1,:))/divx;

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);

            mesh.C = C;

            mesh.dofs = reshape(1:mesh.counts.ndof,2,[])';
            mesh.free = true(mesh.counts.nnod,2);
        end

        function mesh = armar_K(mesh)
            dofselem = 8;
            nnz = mesh.counts.nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
            
            I = zeros(nnz,1);
            J = zeros(nnz,1);
            V = zeros(nnz,1);
            cont = 1;

            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Kel = crearK_Q4(nodos, mesh.C);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

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

            mesh.K = sparse(I,J,V);
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
            draw_Mesh(mesh.elems,mesh.nodes, 'NodeLabel',true,'Type','Q4','Color','b')
            hold off

            % Deformada
            x = mesh.nodes(:,1);
            y = mesh.nodes(:,2);
 
            x_deformada = x + escala*mesh.U(1:dofselem:ndof);
            y_deformada = y + escala*mesh.U(2:dofselem:ndof);
            nodos_deformada = [x_deformada y_deformada];
            
            figure(3)
            draw_Mesh(mesh.elems,mesh.nodes,'Type','Q4','Color','b')
            hold on
            draw_Mesh(mesh.elems,nodos_deformada,'Type','Q4','Color','k')
            hold off
        end
    end


end