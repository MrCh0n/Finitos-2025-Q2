function [] = plot_Q4(nodos, deformaciones, mult)
% nodos son todos los nodos del problema
% deformaciones es la U

plot(nodos(:,1), nodos(:,2), 'k*');
hold on

todo = mult*reshape(deformaciones,2,[])';
deformada = nodos + todo;
plot(deformada(:,1), deformada(:,2), 'bo');
end