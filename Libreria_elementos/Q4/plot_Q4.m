function [] = plot_Q4(nodos,elem, dofs, deformaciones, divisiones, mult)
% nodos son todos los nodos del problema
% deformaciones es la U

nelem = size(elem,1);

plot(nodos(:,1), nodos(:,2), 'k*');
hold on

todo = mult*reshape(deformaciones,2,[])';
deformada = nodos + todo;
plot(deformada(:,1), deformada(:,2), 'bo');
end