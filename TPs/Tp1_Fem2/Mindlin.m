classdef Mindlin < handle
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
        function mesh = Mindlin(L,R,tita,divx,divy,E,v,t)
            dofselem = 6;
            [mesh.nodes, mesh.elems] = mallador_ej1(L/2,R,tita,divx,divy);

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

        function mesh = armar_R(mesh,q)
            dofselem = 6;
            for i=1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);
                dir = mesh.dofs(nodoid,:);
                dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas

                Ae = area(mesh.nodes(nodoid,:));
                
                dir = dir(3:dofselem:4*dofselem); % Carga solamente en z
            
                mesh.R(dir) = mesh.R(dir) + Ae*q/4*ones(4,1); 
            end
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

        function dibujar(mesh)
            dofselem = 6;
            ndof = mesh.counts.ndof;
            % sin deformar
            figure(2)
            draw_Mesh(mesh.elems,mesh.nodes, 'NodeLabel',true,'Type','Q4','Color','b')
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


function Ae = area(nodos)
    v1 = (nodos(2,:)-nodos(1,:))/norm(nodos(2,:)-nodos(1,:));
    
    v3 = cross(v1,(nodos(4,:)-nodos(1,:)));
    
    v3 = v3 / norm(v3); % Normalize the cross product vector
    
    v2 = cross(v3,v1);
    
    Q = [v1;v2;v3]; %v1 es vector fila
    nodosp = (Q*nodos')'; %nodos locales
    nodos = nodosp(:,1:2);

    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(4,1) x1 y1 x1.*y1];
    [w, puntos, n] = gauss([2,2]);
    Ae = 0; % area del elemento
    for i = 1:n
        Neta = [0, 1, 0 puntos(i,2)]/A;
        Nzeta = [0, 0, 1, puntos(i,1)]/A;
        
        D = [Neta; Nzeta];
    
        J = D*nodos;
        mult = abs(det(J))*w(i);
        Ae = Ae + mult;
    end
end