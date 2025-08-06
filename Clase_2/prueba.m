clc;
clear;

valores = struct;

valores.coord = [0 0;
                 6 0;
                 12 0;
                 18 0];

E = 210e9;
I = 700e-6;
valores.elem = [1 2 E I;
                2 3 E I;
                3 4 E 2*I];

valores.BC = [1 2 5];

valores.Fnodos = [1 -30000;
                  2 -30000;
                  3 -60000;
                  5 -30000;
                  6 30000;
                  7 -100000];

vigas(valores);