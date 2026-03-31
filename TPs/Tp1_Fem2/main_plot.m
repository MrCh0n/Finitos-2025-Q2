% plotear

clc
clear all
close all

%% Malla
div = 16; %en que division ver los esfuerzos
divisiones = [2,4,8,16];

wBs = zeros(4,2);
Nxs = zeros(div,2);
Mys = zeros(div,2);
Qys = zeros(div,2);
for i = 1:length(divisiones)
   for type=1:2
         [wB,Nx_p,My_p,Qy_p]= ej1(type,divisiones(i));
         wBs(i, type)=wB;
         if divisiones(i) == div
             Nxs(:,type) = Nx_p;
             Mys(:,type) = My_p;
             Qys(:,type) = Qy_p;
         end
   end
end

figure(1)
hold on; grid on;
x=1:4;
ref_val = -0.3024;
yline(ref_val, '--k', 'LineWidth', 1.5, 'Label', '0,3024', 'LabelHorizontalAlignment', 'left');

plot(x, wBs(:,1), 'o-.b', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');
plot(x, wBs(:,2), 'd-.r', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');


mesh_labels = {'2x2', '4x4', '8x8', '16x16'};
xticks(x);
xticklabels(mesh_labels);
xlabel('Malla');
ylabel('w_B [in]');
title('Convergencia de w_B');

%ylim([-0.55 -0.20]);
set(gca, 'YDir', 'reverse'); % Invertir eje si prefieres ver la deflexión hacia abajo
set(gca, 'Box', 'on', 'TickDir', 'in', 'LineWidth', 1);

% Leyenda
legend('Referencia', 'Cáscara plana Mindlin', 'Cáscara 3D Degenerada', 'Location','northeast');

hold off;

%% Tensiones
tita = 20/div:40/div:40; %los centros de los elementos en el eje tita
tita_ref = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40];
Nx_ref_val = [0, -0.1, -0.3, -0.6, -1.0, -1.4, -1.6, -1.7, -1.7, -1.5, -1.1, -0.4, 0.8, 2.3, 4.0, 5.8, 7.5]*1e4;

figure(2)
hold on; grid on;
plot(tita_ref, Nx_ref_val, 'k', 'LineWidth', 0.5);
plot(tita, Nxs(:,1), 'o-.b', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');
plot(tita, Nxs(:,2), 'd-.r', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');

% Leyenda
xlabel('\theta [º]')
ylabel("N_{x'} [lb/in]")
legend('Referencia', 'Cáscara plana Mindlin', 'Cáscara 3D Degenerada', 'Location','northeast');

hold off;


My_ref_val = [-2080, -2050, -2000, -1920, -1800, -1650, -1480, -1280, -1050, -800, -550, -320, -150, -50, -10, 20, 0];

figure(3)
hold on; grid on;
plot(tita_ref, My_ref_val, 'k', 'LineWidth', 0.5);
plot(tita, Mys(:,1), 'o-.b', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');
plot(tita, Mys(:,2), 'd-.r', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');

% Leyenda
xlabel('\theta [º]')
ylabel("M_{y'} [lb in/in]")
legend('Referencia', 'Cáscara plana Mindlin', 'Cáscara 3D Degenerada', 'Location', 'northeast');

hold off;

Qy_ref_val = [0, 40, 85, 130, 170, 205, 235, 255, 268, 272, 265, 245, 215, 175, 125, 65, 0];

figure(4)
hold on; grid on;
plot(tita_ref, Qy_ref_val, 'k', 'LineWidth', 0.5);
plot(tita, Qys(:,1), 'o-.b', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');
plot(tita, Qys(:,2), 'd-.r', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');

% Leyenda
xlabel('\theta [º]')
ylabel("Q_{y'} [lb/in]")
legend('Referencia', 'Cáscara plana Mindlin', 'Cáscara 3D Degenerada', 'Location', 'northeast');

hold off;
