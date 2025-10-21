%1.3 reticulado, ver guia
clc
clear

sistema = struct('coord',[],'elem',[],'prop',[],'BC',[],'Fnod',[],'Q',[1]);

E = 210e9;
A = 20*20*1e-6;
F = 3000;

sistema.coord = [0 0;
                 2 0;
                 0 2;
                 2 2;
                 4 2];

sistema.elem = [1 2;
                1 3;
                2 3;
                2 4;
                2 5;
                3 4;
                4 5];

for i = 1:size(sistema.elem,1)
    sistema.prop(i,:) = [E A];
end

sistema.BC = [1 5 6];

sistema.Fnod = [8 -F];

barras(sistema);