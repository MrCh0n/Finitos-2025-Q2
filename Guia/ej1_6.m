clc
clear
%TODO hacer esto lindo
 
D = 25; %mm
t = 1.25; %mm
A = (D^2-(D-2*t)^2)*pi/4/1e6; %m^2
I = pi/64*(D^4 -(D-2*t)^4)/1e12; %m^4
E = 200e9; %Pa
b = 0.2; %m
P = 16.4855; %N

mesh.coord = [0;3*b; 3*b; 4*b; 5*b; 5*b; 6*b];

mesh.elem = [1 2;
             3 4;
             4 5;
             6 7];

nnod = size(mesh.coord,1)-1;
ndof_per = 2;
nelem = size(mesh.elem,1);

mesh.dofs = [1:2*2 [3 5] 6:9 [8 10] 11:12];

mesh.dofs = reshape(mesh.dofs,ndof_per,[])';

Libres = true(nnod,ndof_per);
Libres(1,:) = false;
Libres(nnod,:) = false;
Libres = reshape(Libres',1,[]);

R = zeros(ndof_per*nnod,1);
R(6) = -P;

K = zeros(ndof_per*nnod);

for i = 1:nelem
    V = (mesh.coord(mesh.elem(i,2),:) - mesh.coord(mesh.elem(i,1),:));
    L = norm(V);

    Y1 = 12*E*I/L^3;
    Y2 = 6*E*I/L^2;
    Y3 = 4*E*I/L;
    Y4 = 2*E*I/L;

    Kel = [Y1 Y2 -Y1 Y2;
           Y2 Y3 -Y2 Y4;
           -Y1 -Y2 Y1 -Y2;
           Y2 Y4 -Y2 Y3];

    dir = [mesh.dofs(mesh.elem(i,1) , :) mesh.dofs(mesh.elem(i,2) , :)];

    K(dir, dir) = K(dir, dir) + Kel;
end

Rr = R(Libres);

Kr = K(Libres, Libres);

U = zeros(ndof_per*nnod,1);
U(Libres) = Kr\Rr;

U = [U(1:4); U(3); U(5:9); U(8); U(10:12)]