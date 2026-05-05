classdef Q4 < handle
    properties
        nodos
        elems
        bordes
        U
        K
        M
        R
        vibraciones
        material
        campos
        counts
        error
    end

    methods
        function mesh = Q4(bordes,divx,divy,E,v,t,tipo)
            [mesh.nodos.coordenadas, mesh.elems, mesh.bordes] = mallador_cuadrado_Q4(bordes,divx,divy);

            mesh.counts.nnod = size(mesh.nodos.coordenadas,1);
            mesh.counts.ndof = mesh.counts.nnod*2;
            mesh.counts.nelem = size(mesh.elems,1);
            mesh.counts.L = norm(bordes(2,:)-bordes(1,:))/divx;

            mesh.R = zeros(mesh.counts.ndof, 1);
            mesh.U = zeros(mesh.counts.ndof, 1);

            mesh.material.E = E;
            mesh.material.v = v;
            mesh.material.t = t;
            mesh.material.rho = 1;
            
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
            dofselem = 8;
            coord = mesh.nodos.coordenadas;
            elem = mesh.elems;
            dofs = mesh.nodos.dofs;
            C = mesh.material.C;
            funcion = @crearK_Q4;

            mesh.K = armar_K_2D(coord, elem, dofs, dofselem, C, funcion);
        end

        function mesh = armar_R(mesh, i, carga_s, carga_v, type)
                nodoid = mesh.elems(i,:);

                coord = mesh.nodos.coordenadas(nodoid,:);

                Rel = carga_Q4(coord, carga_s, carga_v, "coordenada", type);

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

        function mesh = armar_M(mesh,p,tipo)
            mesh.material.rho = p;
            dofselem = 8;
            coord = mesh.nodos.coordenadas;
            elem = mesh.elems;
            dofs = mesh.nodos.dofs;
            funcion = @masa_Q4;

            mesh.M = armar_M_2D(coord, elem, dofs, dofselem, p, tipo, funcion);
        end

        function mesh = calc_errorzz(mesh)
            deformaciones = zeros(mesh.counts.nelem,3);
            mesh.campos.deformaciones.nodal_bruto = zeros(mesh.counts.nnod, 3);
            count = zeros(mesh.counts.nnod,1);

            for i = 1:mesh.counts.nelem
                nodosid = mesh.elems(i,:);

                coord = mesh.nodos.coordenadas(nodosid,:);

                dir = reshape(mesh.nodos.dofs(nodosid,:)',1,[]);

                Uel = mesh.U(dir);

                epsilon_el = deformaciones_Q4(coord, Uel);

                deformaciones(i,:) = epsilon_el';

                stress_el = stress_elem_a_nodo(coord, epsilon_el');

                mesh.campos.deformaciones.nodal_bruto(nodosid,:) = mesh.campos.deformaciones.nodal_bruto(nodosid,:) + stress_el;
                count(nodosid,:) = count(nodosid,:) + 1;
            end

            if ~isfield(mesh.campos,"Global_matrix")
                mesh.campos.Global_matrix = armar_Global(mesh);
            end

            mesh.campos.deformaciones.bruto = deformaciones;
            mesh.campos.deformaciones.nodal_suavizado = mesh.campos.Global_matrix\mesh.campos.deformaciones.nodal_bruto;
            mesh.campos.deformaciones.nodal_bruto = mesh.campos.deformaciones.nodal_bruto./count;
            
            coord = mesh.nodos.coordenadas;
            elem = mesh.elems;
            bruto = mesh.campos.deformaciones.nodal_bruto;
            suavizado = mesh.campos.deformaciones.nodal_suavizado;
            C = mesh.material.C;

            [zz, U, E] = error_zz_2D(coord, elem, bruto, suavizado, C, @energia_Q4);

            mesh.error.E = E;
            mesh.error.U = U;
            mesh.error.zz = zz;
        end
        
        function mesh = calc_U(mesh)
            free = reshape(mesh.nodos.free', 1, []);
            Kr = mesh.K(free, free);
            Rr = mesh.R(free);

            mesh.U(free) = Kr\Rr;
        end

        function mesh = calc_stress(mesh)
            %TODO poner que use el global devuelta
            coord = mesh.nodos.coordenadas;
            elem = mesh.elems;
            dofs = mesh.nodos.dofs;
            Uel = mesh.U;
            C = mesh.material.C;
            Czz = mesh.material.Czz;
            f_stress = @stress_Q4;
            f_global = @global_Q4;
            f_interpolacion = @elem_a_nodos_Q4;
            [bruto, n_bruto, n_suave] = stress_2D(coord, elem, dofs, Uel, C, Czz, f_stress, f_global, f_interpolacion);

            mesh.campos.stress.bruto = bruto;
            mesh.campos.stress.bruto_nodal = n_bruto;
            mesh.campos.stress.suavizado_nodal = n_suave;
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
            draw_Mesh(mesh.elems,mesh.nodos.coordenadas, 'NodeLabel',true,'Type','Q4','Color','b')
            hold off

            % Deformada
            x = mesh.nodos.coordenadas(:,1);
            y = mesh.nodos.coordenadas(:,2);
 
            x_deformada = x + escala*mesh.U(1:dofselem:ndof);
            y_deformada = y + escala*mesh.U(2:dofselem:ndof);
            nodos_deformada = [x_deformada y_deformada];
            
            figure()
            draw_Mesh(mesh.elems,mesh.nodos.coordenadas,'Type','Q4','Color','b')
            hold on
            draw_Mesh(mesh.elems,nodos_deformada,'Type','Q4','Color','k')
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
                draw_Mesh(mesh.elems,coord,'Type','Q4','Color',color(i))
            end
        end   
    
        function mesh = dibujar_stress(mesh)
            Stress = mesh.campos.stress.bruto;
            %Dibujar graficos
            cmap = jet(256);
            
            ncol = size(cmap,1);
            Stressmax = max(Stress);
            Stressmin = min(Stress);
            
            fig = zeros(7,1);

            for i = 1:7
                fig(i) = figure();
            end

            for i = 1:mesh.counts.nelem
                nodoid = mesh.elems(i,:);%para ordenarlos
            
                Coord = mesh.nodos.coordenadas(nodoid,:);
                
                %sigma xx, yy, xy, zz , I , III vm
                for j=1:7
                    figure(fig(j));
                    dif = (Stressmax(j) - Stressmin(j));
                    if dif > 10
                        t = (Stress(i,j) - Stressmin(j)) / dif;
                    else
                        t = 0.5; % si no hay diferencia de stress
                    end
                    idx = round(1 + t * (ncol-1));
                    plot(polyshape(Coord), 'Facecolor', cmap(idx,:))
                    hold on
                end
                
            end
            e = 1e-16;
            
            
            figure(fig(1))
            colormap(cmap)
            clim([min(Stress(:,1))-e max(Stress(:,1))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{xx} [MPa]');
            hold off
            
            figure(fig(2))
            colormap(cmap)
            clim([min(Stress(:,2))-e max(Stress(:,2))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{yy} [MPa]');
            hold off
            
            figure(fig(3))
            colormap(cmap)
            clim([min(Stress(:,3))-e max(Stress(:,3))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{xy} [MPa]');
            hold off
            
            figure(fig(4))
            colormap(cmap)
            clim([min(Stress(:,4))-e max(Stress(:,4))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{zz} [MPa]');
            hold off
        
            figure(fig(5))
            colormap(cmap)
            clim([min(Stress(:,5))-e max(Stress(:,5))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{I} [MPa]');
            hold off
        
            figure(fig(6))
            colormap(cmap)
            clim([min(Stress(:,6))-e max(Stress(:,6))+e]/1e6);
            colorbar;
            title('Tensión \sigma_{III} [MPa]');
            hold off
            
            figure(fig(7))
            colormap(cmap)
            clim([min(Stress(:,7))-e max(Stress(:,7))+e]/1e6);
            colorbar;
            title("Von Minses");
            hold off
        end
    end
end

function [stress_el] = stress_elem_a_nodo(nodos, stress)
    %% creo el isoparametrico
    cant_puntos = 4;
    
    x1 = [-1; 1; 1; -1];
    y1 = [-1; -1; 1; 1];
    A = [ones(cant_puntos,1) x1 y1 x1.*y1];
    A=inv(A);
    %% Gauss
    [w, puntos, n] = gauss([2,2]);
    
    stress_el=0;
    for i = 1:n
            N = [1, puntos(i,1), puntos(i,2), puntos(i,1)*puntos(i,2)]*A;
    
            Neta = [0, 1, 0 puntos(i,2)]*A;
            Nzeta = [0, 0, 1, puntos(i,1)]*A;
    
            D = [Neta; Nzeta];
    
            J = D*nodos;
      
            mult = w(i);
            Mmin = N'*stress;
            
            stress_el = stress_el + Mmin*mult;
    end% i
end
function M = armar_Global(mesh)
    dofselem = 4;
    nnz = mesh.counts.nelem*dofselem^2;%si ningun nodo se repite se tienen esta cantidad de posibles no zeros
    
    I = zeros(nnz,1);
    J = zeros(nnz,1);
    V = zeros(nnz,1);
    cont = 1;

    for i = 1:mesh.counts.nelem
        nodoid = mesh.elems(i,:);

        coord = mesh.nodos.coordenadas(nodoid,:);

        Mel = global_Q4(coord);

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