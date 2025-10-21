function [nodos, elems, transformada, inicioElem] = mesh_1D(Elems, Nodes, Divisions)
%%funcion para hacer subnodos en un mesh, Elems es la conecciones entre
%%nodos, Nodes son las coordenadas de los nodos y divisions son la
%%cantidad de subelementos que se quieren, minimo 1
%
%Devuelve los nodos nuevos y sus conecciones en nodos y elems
%respectivamente. Transformada es a que indice se van los nodos originales
%y inicio elem son los indices donde empiezan los subelementos que
%pertenecian al mismo elemento

nelem = size(Elems,1);
nnode = size(Nodes,1);

pointsNodeID = zeros(nnode, 1);
Puntocreado = false(nnode,1);

elems = [];
nodos = [];

contador = 1;

inicioElem = zeros(nelem,1);
inicioElem(1) = 1;

transformada = zeros(1,nnode);

for i = 1:nelem
    inicio = Nodes(Elems(i,1), :);
    final = Nodes(Elems(i,2), :);
    
    subnodos = zeros(Divisions(i)+1,3);

    subnodos(:,1) = linspace(inicio(:,1), final(:,1), Divisions(i) +1);
    subnodos(:,2) = linspace(inicio(:,2), final(:,2), Divisions(i) +1);
    subnodos(:,3) = linspace(inicio(:,3), final(:,3), Divisions(i) +1);

    ultimo = contador + Divisions(i);

    subelems = [(contador:ultimo-1)' (contador+1:ultimo)'];

    if Puntocreado(Elems(i,1))
        subnodos(1,:) = [];
        subelems = subelems-1;
        subelems(1,1) = pointsNodeID(Elems(i,1));
    else
        pointsNodeID(Elems(i,1)) = subelems(1,1);
    end

     if Puntocreado(Elems(i,2))
        subnodos(end,:) = [];
        subelems(end,2) = pointsNodeID(Elems(i,2));
    else
        pointsNodeID(Elems(i,2)) = subelems(end,2);
     end
    
    contador = max(max(subelems(:))+1,contador);
    Puntocreado([Elems(i,1) Elems(i,2)]) = true;

    nodos = [nodos; subnodos];
    elems = [elems; subelems];
    
    if ~transformada(Elems(i,1))
        transformada(Elems(i,1)) = max(transformada)+1;
        transformada(Elems(i,1)) = subelems(1,1);
    end

    if ~transformada(Elems(i,2))
        transformada(Elems(i,2)) = max(transformada)+1;
        transformada(Elems(i,2)) = subelems(end,2);
    end

    inicioElem(i+1) = inicioElem(i) + Divisions(i);
end
inicioElem(end) = [];
end