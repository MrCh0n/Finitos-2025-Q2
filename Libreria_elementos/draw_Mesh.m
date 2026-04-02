function draw_Mesh(Elem,Nodos,options)
% Sirve para dibujar el Mesh con las coordenadas en Nodos y Elementos. 
% Para dibujar la deformada se puede reemplazar Nodos con Deformada 


arguments
    Elem (:, :) {mustBeInteger}
    Nodos (:, :) {mustBeNumeric}
    options.Color = 'k';
    options.NodeLabel (1, 1) {mustBeNumericOrLogical} = false;
    options.ElementLabel  (1, 1) {mustBeNumericOrLogical} = false;
    options.FontSize (1, 1) {mustBeInteger} = 18;
    options.MarkerSize (1, 1) {mustBeInteger} = 10;
    options.DisplayName {mustBeText} = '';
    options.LineWidth {mustBeNumeric} = 2;
    options.Type {mustBeText} = '';
end

color = options.Color;
linewidth = options.LineWidth;

[nel,nnel] = size(Elem);

%% Reordenamiento vectorizado
if strcmp(options.Type, 'LST')
    Elem = Elem(:,[1 4 2 5 3 6]);
elseif strcmp(options.Type, 'Q8')
    Elem = Elem(:,[1 5 2 6 3 7 4 8]);
end

hold on

%% LABELS DE NODOS
if options.NodeLabel
    if size(Nodos,2)==2
        text(Nodos(:,1), Nodos(:,2), ...
            string(1:size(Nodos,1)), ...
            'Color',color,'FontSize',options.FontSize);
    else
        text(Nodos(:,1), Nodos(:,2), Nodos(:,3), ...
            string(1:size(Nodos,1)), ...
            'Color',color,'FontSize',options.FontSize);
    end
end

%% DIBUJO DE ELEMENTOS
for k = 1:nel
    
    nodes = Elem(k,:);
    
    coords = Nodos(nodes,:);
    coords = [coords; coords(1,:)]; % cerrar
    
    if size(Nodos,2)==2
        if linewidth == 0
            p = scatter(coords(:,1),coords(:,2), ...
                pi*options.MarkerSize^2/4, ...
                'MarkerFaceColor',color,'MarkerEdgeColor',color);
        else
            p = plot(coords(:,1),coords(:,2), ...
                'LineWidth',linewidth,'Color',color, ...
                'Marker','o','MarkerFaceColor',color, ...
                'MarkerSize',options.MarkerSize);
        end
        
        if options.ElementLabel
            c = mean(coords(1:end-1,:),1);
            text(c(1),c(2),num2str(k),'Color','k');
        end
        
    else
        if linewidth == 0
            p = scatter3(coords(:,1),coords(:,2),coords(:,3), ...
                pi*options.MarkerSize^2/4, ...
                'MarkerFaceColor',color,'MarkerEdgeColor',color);
        else
            p = plot3(coords(:,1),coords(:,2),coords(:,3), ...
                'LineWidth',linewidth,'Color',color, ...
                'Marker','o','MarkerFaceColor',color, ...
                'MarkerSize',options.MarkerSize);
        end
        
        if options.ElementLabel
            c = mean(coords(1:end-1,:),1);
            text(c(1),c(2),c(3),num2str(k),'Color','k');
        end
    end
    
    if ~isempty(options.DisplayName) && k == nel
        set(p,'DisplayName',options.DisplayName)
    else
        set(p,'HandleVisibility','off');
    end
end

grid on
xlabel('x'); ylabel('y');
if size(Nodos,2)==3, zlabel('z'); end
hold off
end