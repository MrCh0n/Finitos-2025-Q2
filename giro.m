function [Q] = giro(punto1, punto2)

x1 = punto2-punto1;

flag1 = 1;
flag2 = 1;

if abs(sum(x1)) == norm(x1)
    if (abs(x1(2)) == 1)
        y1 = x1(2)*[1 0 0];
        x1 = -x1;
    else
        y1 = [0 1 0];
        x1 = sum(x1)*x1;
    end
    flag1 = 0;

    flag2 = 0;
end
if flag1
    y = 1;
    a = 1:3;
    for i = a
        if x1(i) == 0
        b = a;
        b(i) = [];
        A = -x1(b(1))*y/x1(b(2));
        y1(b) = [y A];
        y1(i) = 0;

        flag2 = 0;
        end
    end
end

if flag2
    y = [1 1];
    
    c = -x1(1:2)*y'/x1(3);
    
    y1 = [y c];
end

z1 = cross(x1,y1);

enx = x1/norm(x1);
eny = y1/norm(y1);
enz = z1/norm(z1);

Q = [enx; eny; enz];
end