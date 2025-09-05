function [Q] = giro(punto1, punto2)

x1 = punto2-punto1;
enx = x1/norm(x1);
if(enx(2) == 1)
    y2 = x1 + [-1 0 0];%elige el plano
else
    y2 = x1 + [0 1 0];%elige el plano
end

z1 = cross(x1,y2);%TODO a veces da un - de mas fijarse

y1 = cross(z1,x1);


eny = y1/norm(y1);
enz = z1/norm(z1);

Q = [enx; eny; enz];
end