function [errorzz, U, E] = error_zz_2D(nodos, elems, bruto, suavizado, C, energia)
    nelem = size(elems,1);
    U = 0;
    E = 0;
    for i = 1:nelem
        nodosid = elems(i,:);
    
        coord = nodos(nodosid,:);
    
        e_el = bruto(nodosid,:);
        e2_el = suavizado(nodosid,:)-e_el;
    
        Uel = energia(coord, e_el, C);
    
        Eel = energia(coord, e2_el, C);
        
        U = U + Uel;
        E = E + Eel;
    end
    errorzz = sqrt(E/(E + U));
end