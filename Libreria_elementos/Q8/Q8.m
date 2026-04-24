classdef Q8 < handle
    properties
        nodos
        elems
        bordes
        U
        K
        R
        M
        vibraciones
        material
        campos
        counts
        error
    end

    methods
        function mesh = Q8(bordes,divx,divy,E,v,t,tipo)
            [mesh.nodos.coordenadas, mesh.elems, mesh.bordes] = mallador_cuadrado_Q8(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodos.coordenadas,1);
            mesh.counts.ndof = mesh.counts.nnod*2;
            mesh.counts.nelem = size(mesh.elems,1);
            mesh.counts.L = norm(bordes(2,:)-bordes(1,:))/divx;

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);
    
            mesh.material.E = E;
            mesh.material.v = v;
            mesh.material.t = t;

            if upper(tipo) == "STRESS"
                mesh.material.C = t*E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
                mesh.material.Czz = zeros(1,3);
            else
                mesh.material.C = t*E/((1+v)*(1-2*v))*[1-v v 0;v 1-v 0;0 0 (1-2*v)/2];
                mesh.material.Czz = E*v/(1+v)/(1-2*v)*[1,1,0];
            end

            mesh.nodos.dofs = reshape(1:mesh.counts.ndof,2,[])';
            mesh.nodos.free = true(mesh.counts.nnod,2);
        end

        function mesh = armar_K(mesh)
            dofselem = 16;
            nnz = mesh.counts.nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
            
            I = zeros(nnz,1);
            J = zeros(nnz,1);
            V = zeros(nnz,1);
            cont = 1;

            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);

                coord = mesh.nodos.coordenadas(nodoid,:);

                Kel = crearK_Q8(coord, mesh.material.C);

                dir = reshape(mesh.nodos.dofs(nodoid,:)',1,[]);

                %mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
                for a = 1:dofselem
                    for b = 1:dofselem
                        I(cont) = dir(a);
                        J(cont) = dir(b);
                        V(cont) = Kel(a,b);
                        cont = cont+1;
                    end%b
                end%a
            end%i
            
            mesh.K = sparse(I,J,V); 
        end

        function mesh = armar_R(mesh, i, carga_s, carga_v, type)
                nodoid = mesh.elems(i,:);

                coord = mesh.nodos.coordenadas(nodoid,:);

                Rel = carga_Q8(coord, carga_s, carga_v, "coordenada", type);

                dir = reshape(mesh.nodos.dofs(nodoid,:)',1,[]);

                mesh.R(dir) = mesh.R(dir) + Rel;
        end

        function mesh = cond_borde(mesh, nodos, restricciion)
            if restricciion == 3
                mesh.nodos.free(nodos,:) = false;
            else
                mesh.nodos.free(nodos, restricciion) = false;
            end
        end

        function mesh = armar_M(mesh,p)
            mesh.material.rho = p;
            dofselem = 16;
            nnz = mesh.counts.nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
            
            I = zeros(nnz,1);
            J = zeros(nnz,1);
            V = zeros(nnz,1);
            cont = 1;

            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);

                coord = mesh.nodos.coordenadas(nodoid,:);

                Mel = masa_Q8(coord, mesh.material.rho,1);

                dir = reshape(mesh.nodos.dofs(nodoid,:)',1,[]);

                %mesh.K(dir,dir) = mesh.K(dir,dir) + Kel;
                for a = 1:dofselem
                    for b = 1:dofselem
                        I(cont) = dir(a);
                        J(cont) = dir(b);
                        V(cont) = Mel(a,b);
                        cont = cont+1;
                    end%b
                end%a
            end

            mesh.M = sparse(I,J,V);
        end

        function mesh = calc_U(mesh)
            free = reshape(mesh.nodos.free', 1, []);
            Kr = mesh.K(free, free);
            Rr = mesh.R(free);

            mesh.U(free) = Kr\Rr;
        end

        function mesh = calc_stress(mesh)
            mesh.campos.stress.bruto = zeros(mesh.counts.nelem,7);
            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);
        
                dir = mesh.nodos.dofs(nodoid,:);
                dir = reshape(dir', 1, []); %para que sea un vector leyendo primero columnas
            
                Coord = mesh.nodos.coordenadas(nodoid,:);
                
                Uel = mesh.U(dir);
                
                sigmas = stress_Q8(Coord, Uel, mesh.material.C, mesh.material.Czz);
        
                s1 = sigmas(4);
                s2 = sigmas(5);
                s3 = sigmas(6);
        
                mesh.campos.stress.bruto(i,1:4) = sigmas([1 2 3 6]);
                mesh.campos.stress.bruto(i,5) = max([s1 s2 s3]); %sigma_1
                mesh.campos.stress.bruto(i,6) = min([s1 s2 s3]); %sigma_3
                mesh.campos.stress.bruto(i,7) = sqrt(((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)/2);%von Mises
            end
        end
        
        function mesh = calc_frecuencias(mesh)
            free = reshape(mesh.nodos.free', 1, []);
            ndof = mesh.counts.ndof;
            Kr = mesh.K(free, free);
            Mr = mesh.M(free,free);
            
            k=20;
            [Avec,Aval] = eigs(Kr,Mr,k,'smallestabs');
            nodos_libres = size(Aval,1);
            % Extraer los valores propios y vectores propios
            [f_n,idx]= sort(sqrt(diag(Aval)) / (2 * pi)); % Frecuencias naturales en Hz
            mesh.vibraciones.modos = zeros(nodos_libres,ndof); %menos la cantidad de restricciones
            mesh.vibraciones.modos(:,free) = Avec(:,idx)'; % Modos de vibración
            
            
            mesh.vibraciones.frecuencias = f_n;
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
            figure()
            draw_Mesh(mesh.elems,mesh.nodos.coordenadas, 'NodeLabel',true,'Type','Q8','Color','b')
            hold off

            % Deformada
            x = mesh.nodos.coordenadas(:,1);
            y = mesh.nodos.coordenadas(:,2);
 
            x_deformada = x + escala*mesh.U(1:dofselem:ndof);
            y_deformada = y + escala*mesh.U(2:dofselem:ndof);
            nodos_deformada = [x_deformada y_deformada];
            
            figure()
            draw_Mesh(mesh.elems,mesh.nodos.coordenadas,'Type','Q8','Color','b')
            hold on
            draw_Mesh(mesh.elems,nodos_deformada,'Type','Q8','Color','k')
            hold off
        end
       
        function mesh = dibujar_frecuencias(mesh)
            figure()
            hold on
            color = ['r','g','b','k'];
            for i = [1 2 3]
                x = mesh.vibraciones.modos(i,1:2:end)'+mesh.nodos.coordenadas(:,1);
                y = mesh.vibraciones.modos(i,2:2:end)'+mesh.nodos.coordenadas(:,2);
                coord = [x y];
                draw_Mesh(mesh.elems,coord,'Type','Q8','Color',color(i))
            end
        end
    end


end