function V = seventh (r, R, H, y)
    x=H*r/(R-r);
    ro=(x+y)/x*r;
    V=1/3*pi*(ro.^2.*(x+y)-r^2*x);
end