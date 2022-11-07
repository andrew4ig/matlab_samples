function D = quadroots(a,b,c)
    D= b^2-4*a*c;
    if(D>0)
        fprintf('Уравнение имеет два корня');
        x1=(-b+sqrt(D))/(2*a);
        x2=(-b-sqrt(D))/(2*a);
        fprintf('\nПервый корень:%g',x1);
        fprintf('\nВторой корень:%g',x2);
    end
    if(D==0)
        fprintf('Уравнение имеет один корень');
        x=-b/(2*a);
        fprintf('\nКорень:%g',x);
    end
    if(D<0)
        fprintf('Уравнение не имеет корней');
    end