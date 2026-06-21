function [] = element_check(nodos, elems, check)
    nelem = size(elems,1);
    
    f = figure();
    hold on
    for i = 1:nelem
        nodosid = elems(i,:);

        coord = nodos(nodosid,:);
        
        color = "g";
        
        if ~check(coord)
            color = "r";
        end
        
        plot(polyshape(coord), 'Facecolor', color);
    end
    hold off
end