function [x,y,w] = maxormin (F)
    a=F(1);    b=F(2);    c=F(3);
    x=-b/2/a;
    y=a*x^2+b*x+c;
    w=abs(a)/a;
end