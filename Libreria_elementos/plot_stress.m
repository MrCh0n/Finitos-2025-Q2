function [] = plot_stress(Stress,Elem,Nodos,labels,ind_fig)
%plot_stress grafica las tensiones de los elementos
%   Stress:  nstress x nelem, en cada fila las tensiones de los
%   elementos
%   Elem: elementos del mesh
%   Nodos: Ubicacion de los nodos
%   labels: nombres de las diferentes tensiones
%   ind_fig: a partir de que numero de figura mostrar las tensiones


    nelem = size(Elem,1);
    nstress = size(Stress,2); %cantidad de tensiones a graficar

    %Dibujar graficos
    cmap = jet(256);
    
    ncol = size(cmap,1);
    Stressmax = max(Stress);
    Stressmin = min(Stress);

    for i = 1:nelem
        nodoid = Elem(i,:);%para ordenarlos
    
        Coord = Nodos(nodoid,:);
        
        %sigma xx, yy, xy, zz , I , III vm
        for j=1:nstress
            figure(ind_fig+j)
            dif = (Stressmax(j) - Stressmin(j));
            if dif > 1e-6
                t = (Stress(i,j) - Stressmin(j)) / dif;
            else
                t = 0.5; % si no hay diferencia de stress
            end
            idx = round(1 + t * (ncol-1));
            plot(polyshape(Coord), 'Facecolor', cmap(idx,:))
            hold on
        end
        
    end

    for ifig = 1:nstress
        figure(ind_fig+ifig)
        title(labels{ifig});
        colormap(cmap)
        clim([Stressmin(ifig) Stressmax(ifig)]);
        cb = colorbar;
        xlabel('X Coordinate [m]');
        ylabel('Y Coordinate [m]');
    end

end