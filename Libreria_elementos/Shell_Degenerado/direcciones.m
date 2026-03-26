function [v3_list] = direcciones(nodos,elems)
%calcula el vector v3 promediado de todos los elementos cascara que pertenece
%
%v3 es de 3xnnod

nnod = size(nodos,1);
v3_list = zeros(3,nnod);

for inod = 1:nnod
    [elem_ids,nodo_ids] = find(elems == inod);
    nuso = length(elem_ids);
    v = zeros(nuso,3);
    for iele = 1:nuso
        % para ver que nodo es dentro del elemento (y saber los adyacentes)
        switch nodo_ids(iele)
            case 1
                x1=1; x2=2; y1=1; y2=4;
            case 2
                x1=1; x2=2; y1=2; y2=3;
            case 3
                x1=4; x2=3; y1=2; y2=3;
            case 4
                x1=4; x2=3; y1=1; y2=4;
        end

        %consigo vector perpendicular al elemento iele
        
        nodoselem = nodos(elems(elem_ids(iele),:),:);
        v1 = (nodoselem(x2,:)-nodoselem(x1,:))/norm(nodoselem(x2,:)-nodoselem(x1,:));
    
        v3 = cross(v1,(nodoselem(y2,:)-nodoselem(y1,:)));
    
        v(iele,:) = v3/norm(v3);
    end
    v3 = mean(v,1);
    v3_list(:,inod) = v3/norm(v3);
end
