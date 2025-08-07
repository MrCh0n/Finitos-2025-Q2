%1.1 barra distinto area y material
clc 
clear

sistema = struct('coord',[],'elem',[],'prop',[],'BC',[],'Fnod',[]);

F = 2000;
E1 = 70e9;
E2 = 110e9;
A1 = 50^2*pi/4/10^6;
A2 = A1*30^2/50^2;

sistema.coord = [0 0;
                 2 0;
                 4 0;
                 5 0;
                 8 0]/10;

sistema.elem = [1 2;
                2 3;
                3 4;
                4 5];

sistema.propiedades = [E1 A1;
                        E1 A1;
                        E2 A1;
                        E2 A2];

%TODO porque hay que borrar los Y??
sistema.BC = [1 2 4 6 8 10];

sistema.Fnodos = [3 3*F;
          7 -2*F;
          9 F];

barras(sistema)