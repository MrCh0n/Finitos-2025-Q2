function [w, gp, n] = gauss(ng)
% ng (1, {1 o 2 o 3} ) cantidad de puntos que se quieren por dimension
%
% n cantidad de puntos totales
% w (n,1) pesos para cada punto
% gp (n, {1 o 2 o 3} ) coordenadas de puntos {x o (x,y) o {x,y,z}}

dims = length(ng);
n = prod(ng);
w = ones(n,1);
gp = zeros(n,dims);

switch dims

    case 3
        [w1,gp1] = gauss1D(ng(1));
        [w2,gp2] = gauss1D(ng(2));
        [w3,gp3] = gauss1D(ng(3));

        counter = 1;
        for ig1 = 1:ng(1)
            for ig2 = 1:ng(2)
                for ig3 = 1:ng(3)
                    w(counter) = w(counter)*w1(ig1)*w2(ig2)*w3(ig3);
                    gp(counter,:) = [gp1(ig1) gp2(ig2) gp3(ig3)];
                    counter = counter + 1;
                end
            end
        end

    case 2
        [w1,gp1] = gauss1D(ng(1));
        [w2,gp2] = gauss1D(ng(2));

        counter = 1;
        for ig1 = 1:ng(1)
            for ig2 = 1:ng(2)
                w(counter) = w(counter)*w1(ig1)*w2(ig2);
                gp(counter,:) = [gp1(ig1) gp2(ig2)];
                counter = counter + 1;
            end
        end

    case 1
        [w,gp] = gauss1D(n);
        w = w';
        gp = gp';
                
end



function [w,gp] = gauss1D(n)
%n cantidad de puntos en la cuadratura
%
%w (1,n) pesos de cada punto
%gp (1,n) coordenadas de cada punto
switch n
    case 1
        w  = 2;
        gp = 0;
    case 2
        w  = [1 1];
        a  = sqrt(3)/3;
        gp = [-a a];
    case 3
        w  = [5/9 8/9 5/9];
        a  = sqrt(3/5);
        gp = [-a 0 a];
    case 4
        a  = sqrt((3 - 2*sqrt(6/5))/7);
        b  = sqrt((3 + 2*sqrt(6/5))/7);
        gp = [-b -a a b];
        wa = (18 + sqrt(30))/36;
        wb = (18 - sqrt(30))/36;
        w  = [wb wa wa wb];
    case 5
        a  = 1/3*sqrt(5 - 2*sqrt(10/7));
        b  = 1/3*sqrt(5 + 2*sqrt(10/7));
        gp = [-b -a 0 a b];
        wa = (322 + 13*sqrt(70))/900;
        wb = (322 - 13*sqrt(70))/900;
        w  = [wb wa 128/225 wa wb];
    case 6
        a  = 0.932469514203152;
        b  = 0.661209386466265;
        c  = 0.238619186083197;
        wa = 0.171324492379170;
        wb = 0.360761573048139;
        wc = 0.467913934572691;
        gp = [-a -b -c c b a];
        w  = [wa wb wc wc wb wa];
     case 7
        a  = 0.949107912342759;
        b  = 0.741531185599394;
        c  = 0.405845151377397;
        d  = 0.0;
        wa = 0.129484966168870;
        wb = 0.279705391489277;
        wc = 0.381830050505119;
        wd = 0.417959183673469;
        gp = [-a -b -c d c b a];
        w  = [wa wb wc wd wc wb wa];
    case 8
        a  = 0.183434642495650;
        b  = 0.960289856497536;
        c  = 0.796666477413627;
        d  = 0.525532409916329;
        wa = 0.362683783378362;
        wb = 0.101228536290376;
        wc = 0.222381034453374;
        wd = 0.313706645877887;
        gp = [-a -b -c -d d c b a];
        w  = [wa wb wc wd wd wc wb wa];
     case 9
        a  = 0.324253423403809;
        b  = 0.968160239507626;
        c  = 0.836031107326636;
        d  = 0.613371432700590;
        e  = 0;
        wa = 0.312347077040003;
        wb = 0.081274388361574;
        wc = 0.180648160694857;
        wd = 0.260610696402935;
        we = 0.330239355001260;
        gp = [-a -b -c -d e d c b a];
        w  = [wa wb wc wd we wd wc wb wa];
%     case 10
%         a  = 0.932469514203152;
%         b  = 0.661209386466265;
%         c  = 0.238619186083197;
%         wa = 0.171324492379170;
%         wb = 0.360761573048139;
%         wc = 0.467913934572691;
%         gp = [-a -b -c c b a];
%         w  = [wa wb wc wc wb wa];
end