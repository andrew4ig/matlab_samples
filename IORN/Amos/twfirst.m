function [r, th] = twfirst (r1, th1, r2, th2)
    x1=r1*cosd(th1);
    x2=r2*cosd(th2);
    y1=r1*sind(th1);
    y2=r2*cosd(th2);
    x=x1+x2;
    y=y1+y2;
    r=sqrt(x^2+y^2);
    th=atand(y/x);
end