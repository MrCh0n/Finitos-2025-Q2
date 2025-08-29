function [Q] = giro(punto1, punto2)

x1 = punto2-punto1;

y2 = x1 + [0 1 0];%elige el plano

z1 = cross(x1,y2);%TODO a veces da un - de mas fijarse

y1 = cross(x1,z1);

enx = x1/norm(x1);
eny = y1/norm(y1);
enz = z1/norm(z1);

Q = [enx; eny; enz];
end