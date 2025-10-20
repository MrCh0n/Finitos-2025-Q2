function [] = draw_viga_3D(Nodos, U, mult)
plot3(Nodos(:,1),Nodos(:,2),Nodos(:,3),'ko')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on

mov = reshape(U,6,[])';

Deformada = Nodos + mult*mov(:,1:3);

plot3(Deformada(:,1),Deformada(:,2),Deformada(:,3),'b')
end