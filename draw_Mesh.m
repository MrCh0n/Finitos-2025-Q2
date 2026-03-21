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
    options.Type {mustBeText} = ''; %para decir que es LST o Q8
end

color = options.Color;
linewidth = options.LineWidth;
nel = size(Elem,1);
nnel = size(Elem,2);
xx = zeros(nnel,1);
yy = zeros(nnel,1);
zz = zeros(nnel,1);

%% Cambiar orden con Type
if strcmp(options.Type, 'LST')
    for i=1:nel %cant de elementos
        Elem(i,:) = Elem(i,[1 4 2 5 3 6]);
    end
elseif strcmp(options.Type, 'Q8')
    for i=1:nel %cant de elementos
        Elem(i,:) = Elem(i,[1 5 2 6 3 7 4 8]);
    end
end


hold on


for k = 1:nel
    clear xx
    clear yy
    for i = 1:nnel
        xx(i) = Nodos(Elem(k,i),1);
        yy(i) = Nodos(Elem(k,i),2);
        if options.NodeLabel
            text(xx(i),yy(i),num2str(Elem(k,i)),'VerticalAlignment','bottom','Color',color,'FontSize',options.FontSize);
        end
    end
    xx = [xx xx(1)];
    yy = [yy yy(1)];
    if linewidth == 0
        p = scatter3(xx,yy,zz, pi*options.MarkerSize^2/4, 'Color', color, 'Marker', 'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
    else
        p = plot(xx,yy, 'LineWidth', linewidth,'Color', color, 'Marker', 'o','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',options.MarkerSize);
    end
    if ~isempty(options.DisplayName) && k == nel
        set(p, 'DisplayName', options.DisplayName)
    else
        set(p, 'HandleVisibility', 'off');
    end
    x = mean(Nodos(Elem(k,:),1));
    y = mean(Nodos(Elem(k,:),2));
    if options.ElementLabel
        text(x,y,z,num2str(k),'EdgeColor','k','Color','k','FontSize',18);
    end
end

grid on
% axis equal
xlabel('x')
ylabel('y')
hold off

end