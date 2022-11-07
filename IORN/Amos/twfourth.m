function [th, rad] = twfourth (x,y)
    rad=hypot(x,y);
    th=180*abs(y)/y*(abs(x)/x-1)/(-2)+atand(y/x);
end