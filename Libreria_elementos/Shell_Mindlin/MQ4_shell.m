classdef MQ4_shell < handle
    properties
        nodes
        elems
        dofs
        free
        U
        K
        R
        material
        counts
    end

    methods
        function mesh = MQ4_shell(bordes,divx,divy,E,v,t)
            error("TODO:Falta hacer un mallador")
            %TODO
            dofselem = 6;
            [mesh.nodes, mesh.elems] = mallador_cuadrado_Q4(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodes,1);
            mesh.counts.ndof = mesh.counts.nnod*dofselem;
            mesh.counts.nelem = size(mesh.elems,1);

            mesh.material.E = E;
            mesh.material.v = v;
            mesh.material.t = t;

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);

            mesh.K = spalloc(mesh.counts.ndof, mesh.counts.ndof, mesh.counts.ndof*32);

            mesh.dofs = reshape(1:mesh.counts.ndof,dofselem,[])';
            mesh.free = true(mesh.counts.nnod,dofselem);
        end

        function mesh = armar_K(mesh)
            E = mesh.material.E;
            v = mesh.material.v;
            t = mesh.material.t;
            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);

                nodos = mesh.nodes(nodoid,:);

                Kel = crearK_shellMQ4(nodos, E,v,t);

                dir = reshape(mesh.dofs(nodoid,:)',1,[]);

                mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
            end
        end

        function mesh = armar_R(mesh)
            error("No se tiene un armar R, crearla afuera y usar cargar_R(R)")
            %TODO
            for i=1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);
                dir = mesh.dofs(nodoid,:);
                dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

                %Rel = carga_Shell_MQ4();

                mesh.R(dir) = mesh.R(dir) + Rel; 
            end
        end

        function mesh = cargar_R(mesh,R)
           mesh.R = R;
        end

        function mesh = cond_borde(mesh, nodos, restricciion)
                mesh.free(nodos, restricciion) = false;
        end

        function mesh = calc_U(mesh)
            mesh.free = reshape(mesh.free', 1, []);
            Kr = mesh.K(mesh.free, mesh.free);
            Rr = mesh.R(mesh.free);

            mesh.U(mesh.free) = Kr\Rr;
        end

        function [esfuerzos] = esfuerzo(mesh)
            nelem = mesh.counts.nelem;
            E = mesh.material.E;
            v = mesh.material.v;
            t = mesh.material.t;
            esfuerzos = zeros(nelem,8); %Nx, Ny, Mx, My, Mxy, Qx, Qy
            for i = 1:nelem
                nodoid = mesh.elems(i,:);
            
                dir = mesh.dofs(nodoid,:);
                dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
            
                Coord = mesh.nodes(nodoid,:);

                Uel = mesh.U(dir);
                
                esfuerzos(i,:) = stress_shellMQ4(Coord, Uel, E,v,t);
            end
        end

        function dibujar(mesh)
            dofselem = 6;
            ndof = mesh.counts.ndof;
            % sin deformar
            figure(2)
            %draw_Mesh(mesh.elems,mesh.nodes, 'NodeLabel',true,'Type','Q4','Color','b')
            hold off

            % Deformada
            x = mesh.nodes(:,1);
            y = mesh.nodes(:,2);
            z = mesh.nodes(:,3);
            
            escala = 20;
            x_deformada = x + escala*mesh.U(1:dofselem:ndof);
            y_deformada = y + escala*mesh.U(2:dofselem:ndof);
            z_deformada = z + escala*mesh.U(3:dofselem:ndof);
            nodos_deformada = [x_deformada y_deformada z_deformada];
            
            figure(3)
            draw_Mesh(mesh.elems,mesh.nodes,'Type','Q4','Color','b')
            hold on
            draw_Mesh(mesh.elems,nodos_deformada,'Type','Q4','Color','k')
            hold off
        end
    end


end