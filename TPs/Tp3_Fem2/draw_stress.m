function [] = draw_stress(nodos, elems, bruto)
    Stress = bruto;
    nelem = size(bruto,1);
    %Dibujar graficos
    cmap = jet(256);
    
    ncol = size(cmap,1);
    Stressmax = max(Stress);
    Stressmin = min(Stress);
    
    fig = zeros(7,1);

    for i = 1:7
        fig(i) = figure();
    end

    for i = 1:nelem
        nodoid = elems(i,:);%para ordenarlos
    
        Coord = nodos(nodoid,:);
        
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