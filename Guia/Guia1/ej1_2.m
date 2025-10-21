%1.2 Area variable
clc 
clear

sistema = struct('coord',[],'elem',[],'prop',[],'BC',[],'Fnod',[]);
F = 5000;
E = 70e9;
A0 = 1e-3;
L = 0.5;

nelem = 5;

for i = 0:nelem
    sistema.coord = [sistema.coord; i*L/nelem 0];
end

for i = 1:nelem
    sistema.elem = [sistema.elem; i i+1];
    sistema.prop = [sistema.prop; E A0*(1+(i-0.5)/nelem)];
end

sistema.BC = [2*(1:nelem+1) 2*(nelem+1)-1];

sistema.Fnod = [1 -F];

barras(sistema)